#
# nofoma: Normal forms of matrices
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "JordanNormalform",
Subtitle := "Jordan Normalform of matrices",
Version := "1.0",
Date := "28/11/2025", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Alia",
    LastName := "Bonnet",
    Email := "alia.bonnet@rwth.aachen.de",
    IsAuthor := true,
    IsMaintainer := true,
    Place := "Aachen",
    Institution := "RWTH Aachen"
  )
],

SourceRepository := rec(
  Type := "git",
  URL := "https://github.com/AliaBonnet/JordanNormalform",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := "https://github.com/AliaBonnet/JordanNormalform",
README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),


#SourceRepository := rec(
#    Type := "git",
#    URL := Concatenation( "https://github.com/gap-packages/", ~.PackageName ),
#),
#IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
#PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", ~.PackageName ),
#README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
#PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
#ArchiveURL      := Concatenation( ~.SourceRepository.URL,
#                                 "/releases/download/v", ~.Version,
#                                 "/", ~.PackageName, "-", ~.Version ),
#ArchiveFormats := ".tar.gz .tar.bz2",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  
  "This package computes the Jordan normal form of square matrices over 
  any field that is available in GAP.",

BannerString := Concatenation(
"──────────────────────────────────────────────────────────────────────────\n",
"Loading  JordanNormalform ", ~.Version, ", \n",
"by Alia Bonnet \n",
"──────────────────────────────────────────────────────────────────────────\n"),

PackageDoc := rec(
  BookName  := "JordanNormalform",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Jordan normalform",
),

Dependencies := rec(
  GAP := ">= 4.11",
  NeededOtherPackages := [ ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));

