% cenedese et al. (2013)

epsilon = [0.28 0.37 0.58 0.66 0.69 0.74 0.90 1.04 1.05 1.13 1.13 1.16 1.33 1.35 1.37 1.55 ...
           1.68 1.96 2.07 2.11]';
f = [1.5 2.0 1.5 1.5 2.0 1.5 1.5 1.5 2.0 2.0 1.5 1.5 1.5 1.5 2.0 1.5 1.5 1.5 1.5 1.5]';
gc = [4.2 4.2 2.4 4.5 4.5 2.3 2.4 4.4 4.4 4.4 4.4 2.3 4.2 2.3 4.4 2.2 2.3 2.3 2.3 2.2]'/100;
gv = [0.2
      0.2
      0.4
      1.0
      1.0
      0.3
      1.1
      2.0
      2.0
      3.0
      3.0
      2.1
      4.4
      2.9
      4.4
      4.4
      2.0
      3.1
      4.4
      4.4]/100;

hc = [2.0
      1.9
      2.0
      2.0
      1.9
      1.0
      2.0
      1.8
      1.7
      1.9
      2.0
      2.0
      2.0
      2.0
      1.9
      2.0
      1.0
      1.0
      1.0
      1.0]/100;

hv = [3.7
      4.7
      3.9
      4.0
      4.0
      4.0
      3.5
      4.3
      4.0
      3.5
      3.8
      3.0
      3.4
      2.9
      3.5
      2.4
      3.3
      2.9
      2.3
      2.3]/100;

ToQ = [0.21
       0.27
       0.11
       0.08
       0.06
       0.20
       0.36
       0.23
       0.20
       0.12
       0.08
       0.23
       0.09
       0.57
       0.55
       0.73
       0.36
       0.37
       1.72
       1.24];

doRc = [0, 0, -0.03, 0.56, -0.19, 0.04, -0.30, -0.39, 0.00, 0.59, ...
        -0.08, 0.76, 0.66, 1.16, 0.76, 0.32, 0.20, 0.31, 2.61, ...
        2.67]';
Rc = sqrt(gc .* hc)./f;

d = doRc .* Rc;