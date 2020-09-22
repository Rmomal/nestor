#create package
usethis::create_package("/Users/raphaellemomal/Library/Mobile Documents/com~apple~CloudDocs/Rpkg/VEMtree")
usethis::use_package_doc()
usethis::use_roxygen_md()
devtools::document()
#description and dependences
usethis::edit_file("DESCRIPTION")
usethis::use_gpl3_license("Raphaelle Momal")
usethis::use_package("MASS")
usethis::use_package("parallel")
usethis::use_package("ROCR")
usethis::use_package("dplyr")
usethis::use_package("tidyr")
usethis::use_package("mvtnorm")
usethis::use_package("PLNmodels")
usethis::use_package("sparsepca")
usethis::use_package("huge")
usethis::use_package("Matrix")
usethis::use_package("EMtree")
usethis::use_package("blockmodels")
usethis::use_package("rcdd")
usethis::use_package("ggplot2")
usethis::use_package("gridExtra")
usethis::use_package("magrittr")
usethis::use_package("tibble")
usethis::use_package("reshape2")
usethis::use_package("gmp")
usethis::use_package("stringr")
usethis::use_package("graphics")
usethis::use_package("mclust")
usethis::use_package("useful")
# readme
usethis::use_readme_rmd()
#git link
usethis::use_git()

#coverage
usethis::use_coverage(type="codecov")
#travis
usethis::edit_r_environ()
usethis::use_travis() # inclure [ci skip] dans le message du commit pour pusher sans travis
# pkgdown
usethis::use_pkgdown()
#usethis::edit_r_environ()
usethis::git_vaccinate()
usethis::use_pkgdown_travis()

# ssh deployment key manual:
# use_travis_deploy, but manual see https://gist.github.com/gaborcsardi/68960fb01bec74d89611bb9c2b0b7e5a
key <- openssl::rsa_keygen()
pub_key <- travis:::get_public_key(key)
private_key <- travis:::encode_private_key(key) # for travis key
openssl::write_ssh(pub_key) # for git key


# workflow
devtools::document() # then build and restart
devtools::run_examples()
devtools::check()
usethis::use_vignette("Usage","Usage example")
pkgdown::build_site()

##########
# unitary tests
usethis::use_testthat()
usethis::use_test("test_Infer")
usethis::use_test("Simulate")
usethis::use_test("Initialize")
usethis::use_test("Visualize")
usethis::use_test("ExactComp")
devtools::test()
covr::report()
