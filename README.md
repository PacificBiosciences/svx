<h1 align="center">SVX</h1>

<h3 align="center">Structural variation merging</h3>

SVX is a merging tool for structural variations called from PacBio HiFi sequencing data.

## Early version warning

This is a **very early** pre-release of SVX and is under active development. Expect breaking changes at **all** levels, data formats, and system behavior.

### Limitations

- **Only** the [sawfish](https://github.com/PacificBiosciences/sawfish) SV/CNV calling tool is currently supported
- Support for BNDs is currently very limited
- TR-specific logic is unstable
- CNVs are not yet supported
- Only single sample VCFs are supported

## Availability

- The latest SVX Linux binary is [available here](https://github.com/PacificBiosciences/svx/releases)

## Documentation

- [Basic use](docs/guide.md)
- [Command line interface](docs/cli.md)

## Need help?

If you notice any missing features, bugs, or need assistance with analyzing the
output of SVX, please do not hesitate to [reach out by email](mailto:tmokveld@pacificbiosciences.com)
or open a GitHub issue.

## Support information

SVX is currently in active development and is intended for research use only and not for use
in diagnostic procedures. While efforts have been made to ensure that SVX
lives up to the quality that PacBio strives for, we make no warranty regarding
this software.

As SVX is not covered by any service level agreement or the like, please do
not contact a PacBio Field Applications Scientists or PacBio Customer Service
for assistance with any SVX release. Please report all issues through GitHub
instead. We make no warranty that any such issue will be addressed, to any
extent or within any time frame.

### DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE
PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY
KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES
OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A
PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS
SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO
ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY
REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR
IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
