qme-ng: Quiver Mutation Explorer (Version 2.0.1)
================================

Welcome to qme-ng !

QME is a fast quiver manipulation program intended to:
- find the cardinality of the mutation class of a given quiver
- find the length of the green sequences starting from a given quiver

In order to do so, it uses various optimizations that are specific to the
problem at hand, and which result from a careful study by Grégoire Dupont and
Matthieu Pérotin.

Installation (for macOS)
------------------------------------------------------------------------------
1.If you already have macOS Command Line Tools please skip this step. Otherwise please type "xcode-select --install" and then click "install" in the software popup window.

2.Start "Terminal" from Settings/Utilities and go to the folder where you want to install qme-ng using the cd command. Ex: If you want to install qme-ng in /Users/JohnDoe/Documents then please type "cd /Users/JohnDoe/Documents" and then press enter.

3.Type "./install.sh" and then press enter to complete the installation of all the dependencies (gmp and Boost) as well as compilation & linking of the qme-ng program.

Usage
------------------------------------------------------------------------------
Currently qme-ng is only available using Terminal. In the future we may introduce GUI but this is not at the top of our priority. Here is how to use qme-ng:

1.Start "Terminal" from Settings/Utilities and go to the folder where qme-ng is installed using the cd command. Ex: If qme-ng is in /Users/JohnDoe/Documents then please type "cd /Users/JohnDoe/Documents" and then press enter.

2.Type "./qme-ng" to get instructions.

Note:

1.qme-ng requires an empty line at the end of a quiver/ice quiver file.

2.If you use an ice quiver you should type the entire exchange matrix of the ice quiver, not just the extended exchange matrix of the subquiver containing non-frozen vertices.


Features
------------------------------------------------------------------------------

qme-ng can work in two modes:
  1. Quiver mutation class cardinality
  2. Max green sequences length

qme-ng:
 - uses isomorphisms discrimination in order to speed its exploration;
 - uses a fast isomorphism algorithm ([nauty](http://cs.anu.edu.au/~bdm/nauty/));
 - uses [arbitrary precision arithmetic libraries](http://gmplib.org/), and is thus not limited to your CPU registers size;
 - can read [Bernhard Keller's java application](http://www.math.jussieu.fr/~keller/quivermutation/) files as input;
 - produces exploitable outputs;
 - prints its mutations sequences in a format compatible with [Bernhard Keller's java application](http://www.math.jussieu.fr/~keller/quivermutation/);
 - contains tens of small optimizations to cut as many branches in the exploration tree as fast as possible;
 - is free software (BSD License, see the LICENSE file for more details);
 - is actively maintained (by Ying Zhou)!

Contact
------------------------------------------------------------------------------
Webpage: http://mp-bull.github.com/qme-ng/

Webpage of the forked version: http://mathyingzhou.github.com/qme-ng/

Matthieu Pérotin matthieu.perotin(a)bull.net

Ying Zhou yzhou935(a)brandeis.edu
