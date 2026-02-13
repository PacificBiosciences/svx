//! A thread pool for the rust_htslib::bcf module
//!

// NOTE: There's https://github.com/rust-bio/rust-htslib/blob/master/src/tpool.rs but there's no methods that can use it in the bcf module (bgzf::Reader, bgzf::Writer, bam::Reader, bam::IndexedReader all support it)
//       probably worth extending to cover bcf so we do not have to do this work below:

use crate::utils::util::Result;
use rust_htslib::{bcf, htslib};
use std::ptr::NonNull;

/// Thread pool for parallel BGZF decompression across multiple VCF readers.
/// This allows sharing decompression threads among multiple BCF IndexedReader's, Reader's, and Writer's.
///
/// # Safety
/// 1. Assumes rust-htslib IndexedReader, Reader, Writer layout as of 0.50.x.
/// 2. The pool must outlive all future uses of Reader, Writer, IndexedReader.
///
/// Note: The const assertions below do not fully guarantee layout compatibility; they only detect some incompatible changes in rust-htslib updates.
pub struct HtsThreadPool {
    pool: NonNull<htslib::htsThreadPool>,
}

impl HtsThreadPool {
    pub fn new(n_threads: i32) -> Result<Self> {
        unsafe {
            let mut inner = Box::new(htslib::htsThreadPool {
                pool: std::ptr::null_mut(),
                qsize: n_threads * 2,
            });

            inner.pool = htslib::hts_tpool_init(n_threads);
            if inner.pool.is_null() {
                return Err(crate::svx_error!("Failed to initialize thread pool"));
            }

            Ok(Self {
                pool: NonNull::new(Box::into_raw(inner)).expect("Box::into_raw never returns null"),
            })
        }
    }

    ///
    /// # Safety
    /// 1. Assumes rust-htslib IndexedReader layout as of 0.50.x.
    /// 2. The pool must outlive all future uses of IndexedReader.
    ///
    /// Note: The const assertions below do not fully guarantee layout compatibility; they only detect some incompatible changes in rust-htslib updates.
    pub unsafe fn attach_to_indexed_reader(&self, reader: &mut bcf::IndexedReader) -> Result<()> {
        unsafe {
            #[repr(C)]
            struct IndexedReaderLayout {
                inner: *mut htslib::bcf_srs_t,
                header: std::rc::Rc<bcf::header::HeaderView>,
                current_region: Option<(u32, u64, Option<u64>)>,
            }

            const _: () = {
                assert!(
                    std::mem::size_of::<bcf::IndexedReader>()
                        == std::mem::size_of::<IndexedReaderLayout>()
                );
                assert!(
                    std::mem::align_of::<bcf::IndexedReader>()
                        == std::mem::align_of::<IndexedReaderLayout>()
                );
            };

            let reader_ptr = reader as *mut bcf::IndexedReader as *mut IndexedReaderLayout;

            debug_assert!(
                reader_ptr.is_aligned(),
                "IndexedReader pointer is not properly aligned for IndexedReaderLayout cast"
            );

            let bcf_srs_ptr = (*reader_ptr).inner;
            if bcf_srs_ptr.is_null() {
                return Err(crate::svx_error!("IndexedReader inner pointer is null"));
            }

            if (*bcf_srs_ptr).nreaders != 1 {
                let reader_count = (*bcf_srs_ptr).nreaders;
                return Err(crate::svx_error!(
                    "Expected exactly 1 reader in bcf_srs_t, got {reader_count}"
                ));
            }

            let bcf_sr_ptr = (*bcf_srs_ptr).readers;
            if bcf_sr_ptr.is_null() {
                return Err(crate::svx_error!("Readers array is null"));
            }

            let hts_file_ptr = (*bcf_sr_ptr).file;
            self.attach_to_hts_file(hts_file_ptr, "bcf::IndexedReader")?;
        }

        Ok(())
    }

    ///
    /// # Safety
    /// 1. Assumes rust-htslib Reader layout as of 0.50.x.
    /// 2. The pool must outlive all future uses of Reader.
    ///
    /// Note: The const assertions below do not fully guarantee layout compatibility; they only detect some incompatible changes in rust-htslib updates.
    pub unsafe fn attach_to_reader(&self, reader: &mut bcf::Reader) -> Result<()> {
        #[repr(C)]
        struct ReaderLayout {
            inner: *mut htslib::htsFile,
            header: std::rc::Rc<bcf::header::HeaderView>,
        }

        const _: () = {
            assert!(std::mem::size_of::<bcf::Reader>() == std::mem::size_of::<ReaderLayout>());
            assert!(std::mem::align_of::<bcf::Reader>() == std::mem::align_of::<ReaderLayout>());
        };

        let reader_ptr = reader as *mut bcf::Reader as *mut ReaderLayout;

        debug_assert!(
            reader_ptr.is_aligned(),
            "Reader pointer is not properly aligned for ReaderLayout cast"
        );

        let hts_file_ptr = unsafe { (*reader_ptr).inner };

        unsafe { self.attach_to_hts_file(hts_file_ptr, "bcf::Reader") }
    }

    ///
    /// # Safety
    /// 1. Assumes rust-htslib Writer layout as of 0.50.x.
    /// 2. The pool must outlive all future uses of Writer.
    ///
    /// Note: The const assertions below do not fully guarantee layout compatibility; they only detect some incompatible changes in rust-htslib updates.
    pub unsafe fn attach_to_writer(&self, writer: &mut bcf::Writer) -> Result<()> {
        #[repr(C)]
        struct WriterLayout {
            subset: Option<bcf::header::SampleSubset>,
            header: std::rc::Rc<bcf::header::HeaderView>,
            inner: *mut htslib::htsFile,
        }

        const _: () = {
            assert!(std::mem::size_of::<bcf::Writer>() == std::mem::size_of::<WriterLayout>());
            assert!(std::mem::align_of::<bcf::Writer>() == std::mem::align_of::<WriterLayout>());
        };

        let writer_ptr = writer as *mut bcf::Writer as *mut WriterLayout;

        debug_assert!(
            writer_ptr.is_aligned(),
            "Writer pointer is not properly aligned for WriterLayout cast"
        );

        debug_assert!(
            std::rc::Rc::as_ptr(&unsafe { (*writer_ptr).header.clone() })
                == writer.header() as *const _,
            "Writer header pointer layout mismatch"
        );

        let hts_file_ptr = unsafe { (*writer_ptr).inner };

        unsafe { self.attach_to_hts_file(hts_file_ptr, "bcf::Writer") }
    }

    unsafe fn attach_to_hts_file(
        &self,
        hts_file_ptr: *mut htslib::htsFile,
        target: &str,
    ) -> Result<()> {
        if hts_file_ptr.is_null() {
            return Err(crate::svx_error!("{target} htsFile pointer is null"));
        }

        let ret = unsafe { htslib::hts_set_thread_pool(hts_file_ptr, self.pool.as_ptr()) };
        if ret != 0 {
            return Err(crate::svx_error!(
                "Failed to set thread pool on {target} (error code: {ret})"
            ));
        }

        let pool_ref = unsafe { self.pool.as_ref() };
        let thread_count = if pool_ref.pool.is_null() {
            0
        } else {
            unsafe { htslib::hts_tpool_size(pool_ref.pool) }
        };

        log::trace!("Attached thread pool with {thread_count} threads to {target}");

        Ok(())
    }
}

impl Drop for HtsThreadPool {
    fn drop(&mut self) {
        unsafe {
            let raw = self.pool.as_ptr();
            if !(*raw).pool.is_null() {
                htslib::hts_tpool_destroy((*raw).pool);
            }
            drop(Box::from_raw(raw));
        }
    }
}

unsafe impl Send for HtsThreadPool {}
unsafe impl Sync for HtsThreadPool {}
