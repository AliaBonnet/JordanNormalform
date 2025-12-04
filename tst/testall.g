LoadPackage("JordanNormalform");

TestDirectory(DirectoriesPackageLibrary("JordanNormalform","tst"),
              rec(exitGAP     := true,
                  testOptions := rec(compareFunction := "uptowhitespace",
                                     transformFunction := "removenl") ) );