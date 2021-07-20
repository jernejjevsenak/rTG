Dear CRAN  
After my initial submission, I have received instructions to fix the rTg: 

1. Please omit the redundant "The r package" from the description text. (Please do not start the description with "This package", package name, title or similar.). Please add some more details about the package functionality and implemented methods in your Description text.

The description text is now changed to: 'Methods for comparing different regression algorithms for describing the temporal dynamics of secondary tree growth (xylem and phloem). Users can compare the accuracy of the most common fitting methods usually used to analyse xylem and phloem data, i.e., Gompertz function and General Additive Models (GAMs); and an algorithm newly introduced to the field, i.e., Bayesian Regularised Neural Networks (brnn). The core function of the package is XPSgrowth(), while the results can be interpreted using implemented generic S3 methods, such as plot() and summary().'

2. If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form: authors (year) <doi:...>; authors (year) <arXiv:...>; authors (year, ISBN:...); or if those are not available: <https:...> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

The publication related to our new r package is ready for submission for the journal Ecography. I will add the appropriate reference immediately when the publication is accepted. 

3. '\ dontrun{}' should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in '\ dontrun{}' adds the comment ('# Not run:') as a warning for the user. Does not seem necessary. Please unwrap the examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing.You could also replace '\ dontrun{}' with '\ donttest', if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions.

I have removed '\ dontrun{}' from all examples

4. Otherwise, you can also write some tests.

I have added a few simple tests to check the correct function outputs. However, the package is rather small for now, so when it grows, I will also add more complex tests to ensure the correct functionality.

Thank you for your suggestions.  
Best, Jernej 


##  First submission
* This is a first (fixed) submission of the package rTG (radial Tree Growth).

## Test environments
* local OS X install, R 4.0.4
* rhub R-devel, 32/64 bit (https://builder.r-hub.io/status/original/rTG_0.2.2.tar.gz-2367ca3687e44756845de16183e21fdc)
* rhub Ubuntu Linux 20.04.1 LTS (https://builder.r-hub.io/status/original/rTG_0.2.2.tar.gz-c0823bbbfabb41399489ddfba9a7b318)
* rhub Fedora Linux, R-devel (https://builder.r-hub.io/status/original/rTG_0.2.2.tar.gz-372d7045f7f04ce69b2539ddc8d05de7)
* win-check oldrelease (https://win-builder.r-project.org/a0iHarCH1ygW/00check.log)
* win-check release (https://win-builder.r-project.org/8zN7SoWX1Clb/00check.log)
* win-check devel (https://win-builder.r-project.org/m8nVQtQL876I/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 0 NOTEs

## Downstream dependencies
* Will be checked in all future versions
