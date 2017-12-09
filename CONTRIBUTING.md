# How to contribute
If you are interested in contributing, please give this guide a quick read to make sure the code you write is formatted correctly. 

# Formatting and style
The bulk of the formatting is done with uncrustify, and all you have to do is run the format code target and commit the newly formatted code before you make a
merge request. To do so, if for example you make a change to rdl_dynamics, just run `catkin_make --make-args format_rdl_dynamics`. This will result in all of the
c++ files in the `rdl_dynamics` package to be nicely formatted. Make sure to do that to all packages you modify. There are a few things uncrustify won't catch, so i'll list those here, please abide by them,
 - Classes and structs are camel case with the first letter capitalized, e.g. `MyStruct`. No underscores
 - Functions and class method names are camel cased with the first letter lower cased, e.g. `myClassesMethod()`. No underscores
 - variable names are camel case the same as functions and methods. Lowercase first letter and no underscores
 - Make your variable names and function names as clear as possible, even if it means they are 50 characters long. An exception to this rule would be variables from Featherstone's book. It is ok, and even encouraged, that you use the variable name from the book
 - Const correctness is required as well

Additionally, if you have trouble getting your branch to pass the GitLab pipeline, it is likely due to style errors which are checked 
every build using `cppcheck`. When you run `catkin_make run_tests`, the style of your code will be checked and a results folder will 
be generated in `<path_to_catkin_ws>/build/static_analysis`. That file should tell you if it found any style problems, so if the 
pipeline is failing and you're not sure why you should check there.

# Error handling
- Code must compile with all warnings turned on
- Exceptions are ok. Use the RdlException though and provide a useful message

# Documentation
Any methods, public or private, must be commented with doxygen style comments as well as any member variables created. There are lots of comments in the code, it should be pretty easy to just pattern match something that's already there.

# Bug fixes
If you believe you have found a bug and would like to provide a fix(thanks!), you are required to also submit in your pull request a unit test that exposes the bug you are fixing. If you aren't interested in fixing it yourself, create an issue giving as much information as possible as well as a failure case so the developer(s) can reproduce the bug. Also make sure that your pull request is back onto the develop branch.

# New features
If you are interested in adding a completely new feature, please follow these steps

1. Create an issue so folks know you're interested in the feature and so others can weigh in. This is important because some features the developer(s) may not want in RDL, so it's encouraged that you get the thumbs up to keep you from wasting time. Also others might be able to help or have good implementation ideas.
2. Fork RDL. I think you'll have to put your fork on GitLab for the PR back to upstream to work properly.
3. Create your branch off the `develop` branch. If the issue created for this feature is #123, then please make your feature branch named `feature/RDL-123`.
4. Implement the feature in your fork, unit testing as rigorously as possible all added functionality
5. Once your feature is implemented, tested, and documented and it's ready for prime time, create a merge request to the develop branch and the developer(s) can start reviewing it.
6. Once it has gotten the thumbs up from reviewers and the build has passed in the CI pipeline, it'll get merged in!