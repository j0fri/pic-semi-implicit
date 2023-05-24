me = 1; %electron mass
n0 = 1; %Number of electrons
e = 1; %fundamental charge
Kb = 1; %boltzmann constant

k = 2*pi; %waveparameter
Te = 0.01; %electron temperature

omega_p = sqrt(4*pi*n0*e^2/me);
disp(['Natural frequency: ', num2str(omega_p),', natural period: ', num2str(2*pi/omega_p)]);

v_th = sqrt(Kb*Te/me);
omega = sqrt(omega_p^2+3*k^2*v_th^2);
disp(['Real frequency: ', num2str(omega),', real period: ', num2str(2*pi/omega)]);

debye = sqrt(k*Te/(4*pi*n0*e^2));
realEigen = -0.22*sqrt(pi)*(omega_p/(k*v_th))^3*exp(-1/(2*k^2*debye^2));
disp(['Book Real component of eigenvalue: ', num2str(realEigen)]);

debye = sqrt(k*Te/(4*pi*n0*e^2));
realEigen2 = -pi/2*(omega_p^3/(2*k^2*v_th^2));
disp(['Chatgpt Real component of eigenvalue: ', num2str(realEigen2)]);
