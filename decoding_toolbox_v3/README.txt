========
CONTENTS
========

1. GENERAL
2. HOW TO CITE
3. ADVANTAGES OF THE TOOLBOX
4. FUNCTIONALITY OF THE TOOLBOX AND BASICS
5. INSTALLING THE TOOLBOX
6. COMPILING LIBSVM (GETTING LIBSVM TO WORK)


==========
1. GENERAL
==========

The Decoding Toolbox has been created for classification of structural and 
functional brain images. It is currently optimized for SPM2, SPM5, and 
SPM8, but can also be extended to other brain image analysis tools.

To get started quick, type
> help decoding_example
in your Matlab command window.

If you want a tutorial or create your own decoding script, 
edit the decoding_tutorial.m

If you want a more extended tutorial with many examples, please read our 
publication (see "2. HOW TO CITE").

For more details on all currently available options, type
> help decoding
or
> help decoding_defaults
and for feature_selection
< help decoding_feature_selection

For a more detailed description of the toolbox, see our publication:
http://journal.frontiersin.org/Journal/10.3389/fninf.2014.00088

Please report any bugs to
martin.hebart@bccn-berlin.de or
kai.goergen@bccn-berlin.de

==============
2. HOW TO CITE
==============

If you used the toolbox, the best favor you can do us is to cite it. Please use
the following reference:
Martin N Hebart*, Kai Goergen* and John-Dylan Haynes (2015). The Decoding Toolbox 
(TDT): A versatile software package for multivariate analyses of functional 
imaging data. Front. Neuroinform. 8:88. doi: 10.3389/fninf.2014.00088.
* Martin N Hebart and Kai Goergen contributed equally to this work.

============================
3. ADVANTAGES OF THE TOOLBOX
============================

The advantages of it are:

1. SIMPLICITY: It is very easy to use. Just try out the decoding_example 
on your SPM.mat using leave-one-run out crossvalidation. 
It is one line of code. If you want more detail, work your way through the
decoding_tutorial.m (one page).

2. GENERALITY: It is quite general purpose (you can do searchlight 
decoding, whole brain decoding or ROI decoding with it). In the current 
beta version we implemented SVM classification with libsvm, SV regression 
with libsvm and pattern similarity analysis using voxel pattern 
correlations. It is general purpose in the classifiers that can be used 
and includes feature selection.

3. FLEXIBILITY: It has a well-defined modular structure and can easily be 
set up for all sorts of classification designs. In addition, it can easily 
be extended by your own algorithm (check the HOWTOEXTEND.txt).

4. READABILITY: You should be able to easily read and adjust the code. 
Although sometimes you have to dig for subfunctions if you want to hack 
the toolbox, this structure makes adjustments a lot easier. Everything is 
commented well which should make it easy to find what you want to edit. 
The transparency makes it valuable for programmers who want to adjust the 
code to their needs.

5. SPEED: It is comparably fast (considering that it uses mainly 
uncompiled code) and uses custom-made functions to speed up processing 
(for example running many F-tests in feature selection in matrix format or 
a custom-made correlation function which is up to 20x faster than that
provided by Matlab).

==========================================
4. FUNCTIONALITY OF THE TOOLBOX AND BASICS
==========================================

Typically, in brain image analyses you would like to know whether some 
regional brain activity pattern is significantly activated. In brain 
image classification you are searching for significant information about 
the classified samples. For that reason we optimized the toolbox for users 
who have a group of subjects and would later want to test whether 
significant information is conveyed in the patterns of brain activity. 
As input, you typically have a number of brain images belonging to several 
categories which you would like to classify. As output you get for each 
subject one or several classification volumes (for searchlight analyses) 
or individual values (e.g. mean cross-validation accuracy in ROIs). In a 
second step the statistical significance can be tested easily using your 
brain image analysis toolbox (e.g. second-level analysis in SPM) or simple 
statistics (e.g. a one-sample t-test in Matlab). For simplicity, we set 
chance to 0 as a default and set all other values around 0 (i.e. for 2 
classes and chance level of 50%, values range from -50 to 50).

Important remark: For calculating the mean output value across cross-
validation steps (e.g. runs), we chose to average across all individual
samples across cross-validation steps, rather than first averaging across 
all samples within a cross-validation step and then again averaging across 
all steps. If you would like to weight all decoding steps equally, set 
   cfg.results.setwise= 1;
   and
   cfg.design.set = 1:length(cfg.design.set);
which provides you with a separate output for each cross-validation step.
Then, average over the resulting output images. Alternatively, write your
own routine (see HOWTOEXTEND.txt) 

=========================
5. INSTALLING THE TOOLBOX
=========================

1. ADD IMAGE PROCESSING SOFTWARE (Default: SPM8)
By default, download & "install" SPM8 (i.e. extract it and add it to your 
matlab path). By default, SPM8 is used to read & write brain images.
Alternatively, you can also use SPM2, SPM5, SPM12, or another toolbox, in 
case it is implemented. See HOWTOEXTEND.txt on how to do this.

2. THE TOOLBOX
Simply extract the toolbox add it to your matlab path.

3. RUN & ENJOY :)


============================================
6. COMPILING LIBSVM (GETTING LIBSVM TO WORK)
============================================

Some people experience problems with the mex-files of libsvm that we 
provided. Libsvm uses code that was created in the C programming language 
which can be a lot faster than normal Matlab code. It was translated into 
files that can be read by the computer and Matlab known as mex-files. 
Sometimes, this translation depends on your particular system. In that 
case, you need to re-compile these files yourself.

**************************
Easy fix for Windows users
**************************

If you use windows and the problem is that a module was not found, try 
downloading the missing dll that your mex file depends on here:
http://www.mathworks.de/support/solutions/en/data/1-2RQL4L/index.html
or install a version of Visual C++ that is compatible with your version of 
Matlab.

**************
How to compile
**************

Alternatively, you may want to use your own compilation or a later version 
of libsvm, so you can follow these 7 steps. 

1. Select a compiler:
a) Linux users: get the latest version of gcc.
b) Windows users: get the appropriate compiler for your version of Matlab 
(check the Mathworks website), probably Microsoft Visual C++ 2010 Express.
If Matlab wants an earlier version (e.g. 2008), a later version will work, 
too. Then set the installation path manually to Matlab using
> mex -setup
Would you like mex to locate installed compilers -> n
Select a compiler: -> select the compiler which is most similar to yours.
If a warning message appears and you are asked: "use path anyway?" -> n
Manually enter the pathname which is similar to the suggested one (e.g. if 
the suggestion is "Microsoft Visual Studio 8.0", then use "Microsoft Visual 
Studio 10.0")
c) Mac users: install the gcc compiler by downloading and installing 
Xcode (need to register). If compilation fails, search the internet for 
"Compiling Matlab Mex Files on a Mac".



1. Use the mex-files provided or download libsvm from
http://www.csie.ntu.edu.tw/~cjlin/libsvm/#download

2. Unzip libsvm and navigate to the folder "Matlab".

3. If your version is more recent than 3.12, you can skip this step. 
Otherwise make the following adjustments (on our request the authors of 
libsvm have introduced the quiet mode -q to the svmpredict options)

a) Open svmpredict.c in Matlab and find and delete (or
comment by using // ) the following code:

if(svm_type==NU_SVR || svm_type==EPSILON_SVR)
    {
        mexPrintf("Mean squared error = %g (regression)\n",error/total);
        mexPrintf("Squared correlation coefficient = %g (regression)\n",
            ((total*sumpt-sump*sumt)*(total*sumpt-sump*sumt))/
            ((total*sumpp-sump*sump)*(total*sumtt-sumt*sumt))
            );
    }
    else
  mexPrintf("Accuracy = %g%% (%d/%d) (classification)\n",
      (double)correct/total*100,correct,total);


b) This step is not necessary when option -q is activated, but you should 
do this anyway: Open SVM.cpp in the path below the folder "Matlab" and 
look for

#if 1
static void info(const char *fmt,...)

and replace the 1 by 0.

4. Rename or delete all compiled mex files. Also make sure that no other 
path contains any mex-files with the same name by typing in the Matlab 
command window

> which svmtrain
and
> which svmpredict

There should be no listing.

5. If you use 32bit, then remove all -largeArrayDims from file make.m in 
the folder "Matlab". For both 32bit and 64bit, then just type "make".
If you get only warnings, but no error, then everything went well. 
Run your code using the compiled mex-files and see if everything worked.

6. Troubleshooting: If you still have problems, you may need to try a 
different compiler. If that still doesn't work, then you could in the 
meantime try a different decoding software (e.g. newton_svm).