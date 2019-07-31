# Spherical test
We select a generic nucleus <sup>84</sup>Zr, constrain its
shape to spherical by constraining its quadrupole
moment to zero value
(<code>icstr = 1</code>, with <code>betac = 0.000</code>)
and perform QFAM calculations
for fixed J and variable 0 <= K <= J.


Due to the Wigner-Eckart theorem, spherical nuclei
should exhibit the strength function response invariant to
multipolarity K for fixed J. In this spherical test,
we verify that assertion numerically.


The demonstrated agreement within 7 most significant digits in
strength response function is satisfying, considering all the
numerical integrations being performed. One can obtain even higher
level of agreement if the induced Coulomb interaction is ignored since
it requires numerical integration of highly oscillating function
(Laplacian of the induced proton vector density).



# Note
In this calculation, we set Broyden's iteration tolerance
(variable <code>tol</code> in <code>iter_fam.f</code>) to
a slightly lower value in order to achieve better agreement.
Default value in the original <code>DIRQFAM</code> code is
<code>1.e-5</code>, but in this test we lowered it to
<code>1.e-8</code>.


