# <code>DIRQFAM v2.0.0</code> 

The corresponding <code>DIRQFAM v1.0.0</code>
<a href="https://www.sciencedirect.com/science/article/pii/S0010465520300345">paper</a> [1]
is located in <code>doc</code> directory.

The corresponding <code>DIRQFAM v2.0.0</code> paper [?] will be (is) located in <code>doc</code> directory.

The corresponding <code>DIRHB</code>
<a href="https://www.sciencedirect.com/science/article/abs/pii/S0010465514000836">paper</a> [2]
is located in <code>doc</code> directory.

Tested versions of the <code>DIRQFAM</code> code are located at
<a href="https://github.com/abjelcic/DIRQFAM/releases">releases</a> page. 






# How to use
The <code>DIRQFAM</code> code is built upon the <code>DIRHB</code> program package [2] for
the solution of the stationary relativistic Hartree-Bogoliubov equations for even-even open-shell
nuclei with axially symmetric quadrupole deformation. <code>DIRQFAM</code> complements
<code>DIRHBZ</code> code with the QFAM solver to calculate the multipole response for systems with 
axially symmetric quadrupole deformation. 






## Input
The input data are provided via the <code>dirqfam.dat</code> file, and are separated
into two parts: first part is the same as in Ref. [2] and determines the input parameters
for the ground state calculation, while the second part serves as an interface for QFAM parameters.


<strong>Ground state parameters</strong> (same as in Ref. [2])

* <code>n0f</code>, <code>n0b</code>:<br>
Number of oscillator shells used in expanding the large component of Dirac spinor
(small component is expanded in <code>n0f+1</code> shells) and wave functions of meson fields
respectively. Both <code>n0f</code> and <code>n0b</code> <strong>must be even</strong> numbers. Recommended value of
<code>n0f</code> depends on the nucleus and one should in principle rerun the calculation with
larger <code>n0f</code> and compare the difference in output to establish whether the
convergence is satisfying. Recommended value of <code>n0b</code> is at least <code>2(n0f+1)</code>.
Since one very rarely uses more than <code>n0f=24</code> shells, we recommend fixing the value of <code>n0b=50</code>,
in which case one doesn't have to worry about the <code>n0b</code> parameter.

* <code>beta0</code>, <code>betai</code>:<br>
Deformation parameter of the oscillator basis <code>beta0</code> and of the initial Woods-Saxon potentials 
<code>betai</code> respectively. In order to improve accuracy, these parameters should be close
to the self-consistent ground state quadrupole deformation. For example, when dealing with a nucleus having
deformation parameter β=+0.550, one should use <code>beta0=+0.550</code> and <code>betai=+0.550</code>.

* <code>inin</code>:<br>
The starting parameters for the initial potentials and initial pairing field respectively.
If set to <code>1</code>, the code starts with Woods-Saxon model as initial guess for the
self-consistent potentials and with diagonal pairing field respectively. Otherwise, if set to <code>0</code>,
the code uses the data from previous run stored in <code>dirhb.wel</code> and/or <code>dirhb.del</code>
as starting potentials and pairing field respectively. For casual users, we recommend using the value of <code>1</code>.

* Even-even nuclide to be computed. Element name, followed by the mass number. <strong>If the element name
has only one character, it should begin with an underscore</strong>, eg. <code>\_C 12</code>,
<code>\_O 16</code> and <code>\_U 238</code>. Otherwise, for two character elements simply type e.g. <code>Zr 100</code>.

* <code>Init.Gap</code>:<br>
Initial pairing gap (in MeV) of diagonal pairing field for protons and neutrons, relevant
if <code>inin</code> for pairing field is set to <code>1</code>. We recommend the generic value of 1 MeV.

* <code>Force</code>:<br>
Acronym of the parameter set of the selected energy density functional.
In current version of the code, <code>DD-PC1</code> and <code>DD-ME2</code> are available.

* <code>icstr</code>, <code>betac</code>, <code>cquad</code>:<br>
The quadrupole deformation constraint control parameters. If <code>icstr</code> is set to
<code>0</code>, the quadrupole constraint is not included. If <code>icstr</code> is
set to <code>1</code>, then the constrained value of expected deformation <code>betac</code>
is imposed with the stiffness constant <code>cquad</code>. We recommend the
value of <code>cquad = 0.010</code>, but if the code fails to constrain
the quadrupole moment, one should increase it keeping in mind that too large
stiffness constant may disrupt the convergence of iterations.
The contrained value of deformation <code>betac</code> is defined as:
$$\beta = \sqrt{\frac{5\pi}{9}} \frac{1}{A R_0^2} \int \rho_v(\boldsymbol{r}) (2z^2-r_\perp^2) d\boldsymbol{r} = \frac{4\pi}{3} \frac{1}{A R_0^2} \int \rho_v(\boldsymbol{r}) |\boldsymbol{r}|^2 Y_{20}(\theta,\varphi) d\boldsymbol{r},$$
where $\rho_v(\boldsymbol{r})$ is axially symmetric self-consistent ground state isoscalar-vector density, and $R_0=1.2 A^{1/3}$ fm.

Mind the alignment of input parameters, for example, parameter <code>n0f</code>
should be written in 5 character width after the equality sign. An example of
<code>dirqfam.dat</code> file is provided and the user <strong>must follow the same
alignment pattern</strong>.


<strong>QFAM parameters interface</strong>

* <code>Calculation type</code>:<br>
Value <code>0</code>: Free response is calculated for a given range of energies.
Value <code>1</code>: Self-consistent response is calculated for a given range of energies.
Value <code>2</code>: Self-consistent response is calculated for a given energy and various data are printed.
Value <code>3</code>: Self-consistent solution is calculated along a circular contour and the contour integral is calculated.

* <code>Include Coulomb</code>, <code>Include pairing</code>:<br>
If set to <code>0</code>/<code>1</code>, the Coulomb interaction or pairing is omitted/included
both in ground state and in the QFAM calculation respectively.

* <code>NGH</code>, <code>NGL</code>:<br>
Parameters defining the size of the Gaussian quadrature grid.
<code>NGH</code> is the number of Gauss-Hermite nodes in z>0 direction and
<code>NGL</code> is the number of Gauss-Laguerre nodes in r direction.
One should use at least <code>NGH=max(n0f+1,n0b)</code> and <code>NGL=max(2(n0f+1),n0b)</code>.
We recommend fixing these values to <code>NGH=25</code> and <code>NGL=50</code>, since one rarely uses more than
<code>n0f=24</code> and <code>n0b=50</code> shells.

* <code>Smearing gamma</code>:<br>
The imaginary part of the complex frequency (smearing width) given in MeV.
Reasonable value is around 0.05-0.50 MeV. This parameter is relevant only if calculation
type is se to <code>0</code>, <code>1</code> or <code>2</code>.

* <code>Solver tolerance</code>:<br>
Relative residual error tolerance for GMRES solver. We recommend using the value of <code>1.e-5</code>, which has shown to give
the strength function accurate up to 4 most significant digits.

* <code>Arnoldi vectors</code>:<br>
Maximum number of Arnoldi vectors used in GMRES solver. This is the limit on the number of QFAM solver steps.
We recommend using the value of <code>70</code>. If the GMRES solver fails to satisfy the relative residual error tolerance,
we recommned increasing this value, however keep in mind that this means larger memory consumption of the program. 

* <code>J multipolarity</code>, <code>K multipolarity</code>, <code>Isospin</code>:<br>
Values of J, K, and T, that define the multipole excitation operator. In the current version,
J value is restricted to 0 <= J <= 5. Multipolarity K should be 0 <= K <= J.
<code>Isospin</code> selects isoscalar/isovector excitation if set to
<code>0</code>/<code>1</code> respectively.

* <code>Omega start</code>, <code>Omega end</code>, <code>Delta omega</code>:<br>
Parameters (MeV) that control the starting point, the ending point and
the increment of the energy range over which the response is calculated.
Relevant only if the calculation type is set to <code>0</code> or <code>1</code>.

* <code>Omega print</code>:<br>
The energy for which the self-consistent solution is calculated if calculation type is set to <code>2</code>.

* <code>Omega center</code>, <code>Omega radius</code>, <code>No. points</code>:<br>
Circular contour parameters used if calculation type is set to <code>3</code>.
The contour is a circle centered at <code>Omega center</code> MeV with radius <code>Omega radius</code> MeV.
Number of integration points used for contour integration is selected via <code>No. points</code> parameter.


Mind the input format, all selected values should be aligned between
the equality sign <code>=</code> and the sentinel <code>|</code>.






## Output
The output of the calculation is divided into two parts. The first output file
<code>dirhb.out</code>, located in the <code>GS\_output</code> directory, contains
the information on the ground state calculation. Detailed description of
this file can be found in [2]. 

The second part of the output relevant for the QFAM calculation is located
in the <code>QFAM\_output</code> directory.
The calculated strength function is written in the <code>strength.out</code> file.
If calculation type is set to <code>3</code>, the value of the contour integral is
also printed in <code>strength.out</code>.
If calculation type is set to <code>2</code>, <code>QFAM\_output</code> contains additional details
and information about the self-consistent solution obtained for a given energy <code>Omega print</code>.






## How to run
The programming language of the <code>DIRQFAM</code> code is Fortran and the user
should provide an implementation of the
<a href="http://www.netlib.org/blas/">BLAS</a> and
<a href="http://www.netlib.org/lapack/">LAPACK</a> (version 3.6.0. or higher)
linear algebra libraries.
Since the code depends heavily on <code>zgemm</code>, <code>dgemm</code> and
<code>dgemv</code> subroutines, the user should provide an efficient implementation
of the BLAS library. We recommend an open source implementation
<a href="https://www.openblas.net/">OpenBLAS</a>, or freely available
<a href="https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html">Intel<sup>®</sup> oneAPI Math Kernel Library</a>
as a part of the Intel<sup>®</sup> oneAPI Base Toolkit.

The code is compiled by standard Makefile build automation which is set to work with
the <a href="https://gcc.gnu.org/fortran/">GFortran</a> compiler.

If the user invokes <code>make</code> command, the 
compilation of the code will produce the executable file <code>run</code>.
The code is then executed by invoking the <code>./run</code> command.
If the user invokes <code>make dbg</code>, the code is compiled in debug mode
in which various additional checks are performed, and the executable file
<code>dbg</code> is generated. Since the executable produced
in debug mode is considerably slower, this mode should be used only for
testing and developing.
If OpenBLAS is employed, the command <code>export OPENBLAS\_NUM\_THREADS=4</code> can be invoked
to select the number of threads used by OpenBLAS.
If Intel<sup>®</sup> oneAPI Math Kernel Library is employed, the command <code>export MKL\_NUM\_THREADS=4</code> can be invoked
to select the number of threads used by the Intel<sup>®</sup> oneAPI Math Kernel Library.





## Test examples
In <code>test</code> directory we provide a set of input files and expected output of the program.
Each test contains an <code>input</code> and <code>output</code> subdirectory which contain the input and expected output data.
It is advisable that the user first tries to reproduce these values.






# Performance benchmark
Benchmark was done using Intel<sup>®</sup> Core<sup>®</sup> i7-9750H @ 2.60GHz machine (laptop).
BLAS and LAPACK are provided via Intel<sup>®</sup> oneAPI Math Kernel Library, which are forced to
run using a single thread, i.e. the entire benchmark is performed using a single thread.

We select the isoscalar J=5, K=3 excitation with
Gaussian quadrature grid: <code>NGH=25</code>, <code>NGL=50</code>.
The <code>DD-ME2</code> parametrization is used with <code>n0b=50</code> shells, and 70
Arnoldi vectors stored in the memory are used by the GMRES solver.

The following table shows running time per QFAM iteration together
with total memory usage. It takes roughly 30-60 iterations to reach self-consistency for a given excitation
energy, depending on the self-consistency tolerance.

| <code>n0f</code> | Memory[GB] | Time[s] | 
| :--------------: | :--------: | :-----: | 
| 10               | 0.66       | 0.25    |
| 12               | 0.90       | 0.52    |
| 14               | 1.31       | 1.04    |
| 16               | 1.94       | 1.98    |
| 18               | 2.99       | 3.65    |
| 20               | 4.54       | 6.48    |






# References
[1] A. Bjelčić, T.Nikšić, Comp. Phys. Comm. 253, 107184 (2020).

[2] T. Nikšić, N. Paar, D. Vretenar, P. Ring, Comp. Phys. Comm. 185, 1808 (2014).
