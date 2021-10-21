setwd("Z:\\2021\\Javan\\candidaPackage")
# There isn't an OrgDb for C. auris on the AnnotationHub, so you will need to generate one yourself, using makeOrgPackageFromNCBI in the AnnotationForge package. You can emulate the example in ?makeOrgPackageFromNCBI, substituting Candida and auris for genus and species, and 498019 for the tax_id. You could also use your actual email and name, but it doesn't really matter so long as the maintainer field has a name and email, and the email is bracketed by < and >. If you don't do that the package won't install.
# Make sure you have the current version of AnnotationForge (AnnotationForge_1.35.2). I patched a bug recently that would cause it to error out on a species like this. And be prepared to wait for a while - that function downloads and parses a huge amount of data.
# Once you have the package built you can then do install.packages("org.Cauris.eg.db", repos = NULL) and if you are on Windows, add a type = "source" to the install.packages arguments

#BiocManager::install(c("AnnotationHub","AnnotationForge"))

library(AnnotationHub)
library(AnnotationForge)

#Make the organism db from NCBI
makeOrgPackageFromNCBI(version = "1.0.0",
                       author = "Javan Okendo <javanokendo@gmail.com>",
                       maintainer = "Javan Okendo <javanokendo@gmail.com>",
                       outputDir = ".",
                       tax_id = "498019",
                       genus = "Candida",
                       species = "auris",
                       NCBIFilesDir=".")
