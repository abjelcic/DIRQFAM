# DIRQFAM
<i>Implementation of the quasiparticle finite amplitude method
within the relativistic self-consistent mean-field framework:
the program DIRQFAM</i>.

The corresponding <code>DIRQFAM</code> paper [1]  (DD-PC1) can
be found at the following <a href="https://arxiv.org/">link (arXiv)</a>.</br>
The corresponding <code>DIRQFAM</code> paper [2]  (DD-ME2) can
be found at the following <a href="https://arxiv.org/">link (arXiv)</a>.</br>




# How to use
The <code>DIRQFAM</code> code is built upon the <code>DIRHB</code> program package [3] for
the solution of the stationary relativistic Hartree-Bogoliubov equations for even-even open-shell
nuclei with axially symmetric quadrupole deformation. <code>DIRQFAM</code> complements
<code>DIRHBZ</code> code with the QFAM solver to calculate the multipole response for systems with 
axially symmetric quadrupole deformation. 





## Input
The input data are provided in <code>dirqfam.dat</code> file, and are separated
into two parts: first part is the same as in [3] and determines the input parameters
for ground state calculation, while the second part serves as an interface for QFAM parameters.



<strong>Ground state parameters</strong> (same as in [3])


* <code>n0f</code>, <code>n0b</code>:</br>
Number of oscillator shells used in expanding the large component of Dirac spinor
(small component is expanded in <code>n0f+1</code> shells) and wave functions of meson fields
respectively. Both <code>n0f</code> and <code>n0b</code> must be even number. Recommended value of
<code>n0f</code> depends on the nucleus and one should in principle rerun the calculation with
larger <code>n0f</code> and compare the difference in output to establish whether the
convergence is satisfying. Recommended value of <code>n0b</code> is at least <code>2(n0f+1)</code>.


* <code>beta0</code>, <code>betai</code>:</br>
Deformation parameter of the oscillator basis and of the initial Woods-Saxon potentials
respectively. In order to improve accuracy, these parameters should be close
to the self-consistent ground state quadrupole deformation.


* <code>inin</code>:</br>
The starting parameters for the initial potentials and initial pairing field respectively.
If set to <code>1</code>, the code starts with Woods-Saxon as initial guess for the
potentials and with diagonal pairing field respectively. Otherwise, if set to <code>0</code>,
the code uses the data from previous run stored in <code>dirhb.wel</code> and/or <code>dirhb.del</code>
as starting potentials and pairing field respectively.


* Even-even nuclide to be computed. Element name, followed by the mass number. If the element name
has only one character, it should begin with an underscore, eg. <code>_C 12</code>,
<code>_O 16</code> and <code>_U 238</code>.


* <code>Init.Gap</code>:</br>
Initial pairing gap (in MeV) of diagonal pairing field for protons and neutrons, relevant
if <code>inin</code> for pairing field is set to <code>1</code>. The default value is 1 MeV.


* <code>Force</code>:</br>
Acronym of the parameter set of the selected energy density functional.
In current version of the code, <code>DD-PC1</code> and <code>DD-ME2</code> are available.


* <code>icstr</code>, <code>betac</code>, <code>cquad</code>:</br>
The quadrupole constraint control parameters. If <code>icstr</code> is set to
<code>0</code>, the quadrupole constraint is not included. If <code>icstr</code> is
set to <code>1</code>, then constrained value of expected deformation <code>betac</code>
is imposed with the stiffness constant <code>cquad</code>. We recommend the
default value of <code>cquad = 0.100</code>, but if the code fails to constrain
the quadrupole moment, one should increase it keeping in mind that too large
stiffness constant may disrupt the convergence of iterations.


Mind the alignment of input parameters, for example, parameter <code>n0f</code>
should be written in 5 character width after the equality sign. An example of
<code>dirqfam.dat</code> file is provided and the user should follow the same
alignment pattern.

</br>


<strong>QFAM parameters interface</strong>


* <code>Calculation type</code>:</br>
Value <code>0</code>: Free response is calculated for a given range of energies.</br>
Value <code>1</code>: Self-consistent response is calculated for a given range of energies.</br>
Value <code>2</code>: Self-consistent response is calculated for a given energy and the
induced density is printed.


* <code>Include Coulomb</code>, <code>Include pairing</code>:</br>
If set to <code>0</code>/<code>1</code>, Coulomb interaction or pairing is omitted/included
both in ground state and the QFAM calculation respectively.


* <code>NGH</code>, <code>NGL</code>:</br>
Number of Gauss-Hermite/Gauss-Laguerre nodes in z>0/r direction respectively.
Recommended values are at least <code>2(n0f+1)</code>.


* <code>J multipolarity</code>, <code>K multipolarity</code>, <code>Isospin</code>:</br>
Values of J, K, and T, that define the multipole operator. In the current version,
J value is restricted to 0 <= J <= 3. Multipolarity K should be 0 <= K <= J.
<code>Isospin</code> selects isoscalar/isovector excitation if set to
<code>0</code>/<code>1</code> respectively.


* <code>Gamma smear</code>:</br>
The smearing width (in MeV) used in the QFAM calculation. Reasonable value is around 0.25-1.00 MeV.


* <code>Omega start</code>, <code>Omega end</code>, <code>Delta omega</code>:</br>
Parameters (in MeV/hbar) that control the starting point, the ending point and
the increment of the energy range over which the response is calculated.
Relevant only if the calculation type flag is set to <code>0</code> or <code>1</code>.


* <code>Omega print</code>:</br>
If the user chooses to calculate the response and print the induced density
for some particular value of energy (calculation type flag set to <code>2</code>),
this value of energy is also provided in the input file.


Mind the input format, all selected values should be aligned between
the equality sign <code>=</code> and the sentinel <code>|</code>.





## Output
The output of the calculation is divided into two parts. The first output file
<code>dirhb.out</code> located in the <code>GS\_output</code> directory contains
the information on the ground state calculation. Detailed description of
this file can be found in [3]. 


The second part of the output relevant for the QFAM calculation is located
in the <code>QFAM\_output</code> directory. The calculated strength function is
written to the <code>strength.out</code> file. If the calculation type flag is
set to 2, an additional file <code>rhov.out</code> is generated which contains
the induced and ground state (vector) density for selected energy. The values printed in the 
<code>rhov.out</code> file are suitable for visualizing the oscillation of time
dependent density. An example of animation easily obtainable from <code>rhov.out</code>
file can be found <a href="http://web.studenti.math.pmf.unizg.hr/~abjelcic/anim.zip">here</a>.
An example of response function easily obtainable from <code>strength.out</code> file can
be found <a href="http://web.studenti.math.pmf.unizg.hr/~abjelcic/anim.zip">here</a>.


## How to run
The programming language of the <code>DIRQFAM</code> code is Fortran and the user
should provide an implementation of the
<a href="http://www.netlib.org/blas/">BLAS</a> and
<a href="http://www.netlib.org/lapack/">LAPACK</a> (version 3.6.0. or higher)
linear algebra libraries.
Since the code depends heavily on <code>zgemm</code>, <code>dgemm</code> and
<code>dgemv</code> subroutines, the user should provide an efficient implementation
of the BLAS library. We recommend an open source implementation
<a href="https://www.openblas.net/">OpenBLAS</a>.
</br>

The code is compiled by standard Makefile build automation which is set to work with
the <a href="https://gcc.gnu.org/fortran/">GFortran</a> compiler.
</br>

If the user invokes <code>make run</code>, first an auxiliary code
will generate <code>dirqfam.par</code> file that contains the relevant information
about the dimensions of various arrays used in the code followed by the
compilation of the code which produces the executable file <code>run</code>.
The code is then executed by invoking the <code>./run</code> command.
</br>
Since the dimensions of various arrays used are statically precalculated, the code
won't rerun if some of the parameters are changed (between two runs) and requires recompilation. This
cumbersome way of using the code will be removed in next upgrades since
the <code>DIRQFAM</code> code is built upon the legacy code <code>DIRHB</code>
written in FORTRAN77 (lack of dynamic allocation).




# Running time benchmark
Benchmark was done on Intel<sup>®</sup> Xeon<sup>®</sup> CPU X5650 @ 2.67 GHz.

Test was performed for J = 3, K = 1 excitation of <sup>20</sup>Ne nucleus with
very dense Gaussian mesh: <code>NGH = 48</code>, <code>NGL = 48.</code><br />

In fact, running time per iteration and memory used depend only on the
multipolarity K, number of nodes in Gaussian mesh, number of shells N in
expansion and number of vectors retained in Broyden's mixing procedure
(default value is 20).


The following table shows running time of subroutines for every QFAM iteration.
It takes roughly 30-60 iterations for a given excitation energy and it depends on
many factors, some of which are the nucleus under consideration, smearing parameter,
vicinity of resonance, density of Gaussian mesh (inaccurate numerical integration
introduces noise and may perturb convergence of iterations) and tolerance of Broyden's
iterations. Therefore, within 5-10 minutes, one can get a 14 shells response for
selected energy. If properly parallelized, provided enough memory and independent
nodes, within few hours one can obtain full 20 shells response of extremely heavy
deformed nucleus with very dense Gaussian mesh.

| N shells  | Total time per iteration | Memory used | fam_drhodkappa | fam_ddensdcurr | fam_dpotentials | fam_dh1 | fam_ddelta | fam_broyden |
| :-------: | :----------------------: | :--------:  | :------------: | :------------: | :-------------: | :-----: | :--------: | :---------: |
| 8         | 0.86 s                   | 0.5 GiB     | 0.12 s         | 0.41 s         | 0.05 s          | 0.24 s  | 0.01 s     | 0.03 s      | 
| 10        | 1.87 s                   | 0.8 GiB     | 0.23 s         | 0.88 s         | 0.05 s          | 0.60 s  | 0.03 s     | 0.08 s      | 
| 12        | 3.79 s                   | 1.5 GiB     | 0.45 s         | 1.80 s         | 0.05 s          | 1.24 s  | 0.07 s     | 0.18 s      | 
| 14        | 7.38 s                   | 2.5 GiB     | 0.82 s         | 3.54 s         | 0.05 s          | 2.37 s  | 0.14 s     | 0.46 s      | 
| 16        | 14.02 s                  | 4.3 GiB     | 1.57 s         | 7.10 s         | 0.05 s          | 4.22 s  | 0.27 s     | 0.81 s      | 
| 18        | 24.49 s                  | 7.1 GiB     | 2.70 s         | 12.57 s        | 0.05 s          | 7.11 s  | 0.60 s     | 1.46 s      | 
| 20        | 41.03 s                  |11.5 GiB     | 4.75 s         | 21.67 s        | 0.05 s          | 11.37 s | 1.12 s     | 2.07 s      | 

<br>
The following table depicts the dependence of running time per QFAM iteration on
multipolarity K for fixed J = 3. Again, test was performed on <sup>20</sup>Ne
nucleus with very dense Gaussian mesh: <code>NGH = 48</code> and <code>NGL = 48</code>.

| N shells | Total time per iteration (K=0) | Total time per iteration (K=1) | Total time per iteration (K=2) | Total time per iteration (K=3) |
| :------: | :----------------------------: | :----------------------------: | :----------------------------: | :----------------------------: |
| 8        | 0.58 s                         | 0.86 s                         | 0.82 s                         | 0.77 s                         |
| 10       | 1.20 s                         | 1.87 s                         | 1.81 s                         | 1.73 s                         |
| 12       | 2.38 s                         | 3.79 s                         | 3.68 s                         | 3.64 s                         | 
| 14       | 4.52 s                         | 7.38 s                         | 7.24 s                         | 7.03 s                         |
| 16       | 8.40 s                         | 14.0 s                         | 13.9 s                         | 13.6 s                         |
| 18       | 14.8 s                         | 24.5 s                         | 24.4 s                         | 23.8 s                         |
| 20       | 25.4 s                         | 41.0 s                         | 41.0 s                         | 41.0 s                         |





<br><br>
To conceive a sense of numerical complexity we are dealing with, the
following table displays the number of vectors (only simplex +i part)
of simplex-y deformed harmonic oscillator basis as a function of number
of shells. It represents typical size of a matrix appearing in QFAM equations.

| N shells | Large <i>f</i> component (N+1)(N+2)(N+3)/6 | Small <i>g</i> component (N+2)(N+3)(N+4)/6 | Total size |
| :------: | :----------------------------------------: | :----------------------------------------: | :--------: |
| 6        | 84                                         | 120                                        | 204        |
| 8        | 165                                        | 220                                        | 385        |
| 10       | 286                                        | 364                                        | 650        |
| 12       | 455                                        | 560                                        | 1015       |
| 14       | 680                                        | 816                                        | 1496       |
| 16       | 969                                        | 1140                                       | 2109       |
| 18       | 1330                                       | 1540                                       | 2870       |
| 20       | 1771                                       | 2024                                       | 3795       |



# References
[1] A. Bjelčić, T.Nikšić, <i>Implementation of the quasiparticle finite
amplitude method within the relativistic self-consistent mean-field framework:
the program DIRQFAM</i>, arXiv preprint submitted to Computer Physics Communications

[2] T. Nikšić, N. Paar, D. Vretenar, P. Ring, <i>DIRHB - a relativistic
self-consistent mean-field framework for atomic nuclei</i>, Comp. Phys.
Comm. 185, 1808 (2014).



