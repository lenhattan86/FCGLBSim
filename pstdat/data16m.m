% data16m.m (16 Machine System Data)
% This is a 16-machine system with 86 transmission lines and
% 68 buses. Data are extracted from the GE final report
% entitled "Singular Perturbations, Coherency and
% Aggregation of Dynamic Systems," pp.6-42, July 1981.
% detailed generator models
% dc exciters on generators 1 to 9
% constant current active loads

% Bus data format
% bus: 
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu),
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%       bus_type - 1, swing bus
%                - 2, generator bus (PV bus)
%                - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)


bus = [...
   1 1.00    0.00   0.00   0.00   2.527  1.1856  0.00 0.00  3  0  0;
   2 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   3 1.00    0.00   0.00   0.00   3.22   0.02    0.00 0.00  3  0  0;
   4 1.00    0.00   0.00   0.00   5.00   1.840   0.00 0.00  3  0  0;
   5 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   6 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   7 1.00    0.00   0.00   0.00   2.34   0.84    0.00 0.00  3  0  0;
   8 1.00    0.00   0.00   0.00   5.22   1.77    0.00 0.00  3  0  0;
   9 1.00    0.00   0.00   0.00   1.04   1.25    0.00 0.00  3  0  0;
   10 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   11 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   12 1.00    0.00   0.00   0.00   0.09   0.88    0.00 0.00  3  0  0;
   13 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   14 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   15 1.00    0.00   0.00   0.00   3.200  1.5300  0.00 0.00  3  0  0;
   16 1.00    0.00   0.00   0.00   3.290  0.32    0.00 0.00  3  0  0;
   17 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   18 1.00    0.00   0.00   0.00   1.58   0.30    0.00 0.00  3  0  0;
   19 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   20 1.00    0.00   0.00   0.00   6.800  1.03    0.00 0.00  3  0  0;
   21 1.00    0.00   0.00   0.00   2.740  1.15    0.00 0.00  3  0  0;
   22 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   23 1.00    0.00   0.00   0.00   2.480  0.85    0.00 0.00  3  0  0;
   24 1.00    0.00   0.00   0.00   3.09  -0.92    0.00 0.00  3  0  0;
   25 1.00    0.00   0.00   0.00   2.24   0.47    0.00 0.00  3  0  0;
   26 1.00    0.00   0.00   0.00   1.39   0.17    0.00 0.00  3  0  0;
   27 1.00    0.00   0.00   0.00   2.810  0.76    0.00 0.00  3  0  0;
   28 1.00    0.00   0.00   0.00   2.060  0.28    0.00 0.00  3  0  0;
   29 1.00    0.00   0.00   0.00   2.840  0.27    0.00 0.00  3  0  0;
   30 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   31 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   32 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   33 1.00    0.00   0.00   0.00   1.12   0.00    0.00 0.00  3  0  0;
   34 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   35 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   36 1.00    0.00   0.00   0.00   1.02  -0.1946  0.00 0.00  3  0  0;
   37 1.00    0.00   0.00   0.00  60.00   3.00    0.00 0.00  3  0  0;
   38 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   39 1.00    0.00   0.00   0.00   2.67   0.126   0.00 0.00  3  0  0;
   40 1.00    0.00   0.00   0.00   0.6563 0.2353  0.00 0.00  3  0  0;
   41 1.00    0.00   0.00   0.00  10.00   2.50    0.00 0.00  3  0  0;
   42 1.00    0.00   0.00   0.00  11.50   2.50    0.00 0.00  3  0  0;
   43 1.00    0.00   0.00   0.00   0.00   0.00    0.00 0.00  3  0  0;
   44 1.00    0.00   0.00   0.00   2.6755 0.0484  0.00 0.00  3  0  0;
   45 1.00    0.00   0.00   0.00   2.08   0.21    0.00 0.00  3  0  0;
   46 1.00    0.00   0.00   0.00   1.507  0.285   0.00 0.00  3  0  0;
   47 1.00    0.00   0.00   0.00   2.0312 0.3259  0.00 0.00  3  0  0;
   48 1.00    0.00   0.00   0.00   2.4120 0.022   0.00 0.00  3  0  0;
   49 1.00    0.00   0.00   0.00   1.6400 0.29    0.00 0.00  3  0  0;
   50 1.00    0.00   0.00   0.00   1.00  -1.47    0.00 0.00  3  0  0;
   51 1.00    0.00   0.00   0.00   3.37  -1.22    0.00 0.00  3  0  0;
   52 1.00    0.00   0.00   0.00  24.70   1.23    0.00 0.00  3  0  0;
   53 1.045   0.00   2.50   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   54 0.98    0.00   5.45   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   55 0.983   0.00   6.50   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   56 0.997   0.00   6.32   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   57 1.011   0.00   5.052  0.00   0.00   0.00    0.00 0.00  2  999  -999;
   58 1.050   0.00   7.00   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   59 1.063   0.00   5.60   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   60 1.03    0.00   5.40   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   61 1.025   0.00   8.00   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   62 1.010   0.00   5.00   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   63 1.000   0.00  10.000  0.00   0.00   0.00    0.00 0.00  2  999  -999;
   64 1.0156  0.00  13.50   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   65 1.011   0.00  35.91   0.00   0.00   0.00    0.00 0.00  1  0  0;
   66 1.00    0.00  17.85   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   67 1.000   0.00  10.00   0.00   0.00   0.00    0.00 0.00  2  999  -999;
   68 1.000   0.00  40.00   0.00   0.00   0.00    0.00 0.00  2  999  -999];


% Line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio, phase shift(deg)

line = [...
   1   2  0.0035  0.0411  0.6987  0      0.;
   1  30  0.0008  0.0074  0.48    0      0.;
   2   3  0.0013  0.0151  0.2572  0      0.;
   2  25  0.007   0.0086  0.146   0      0.;
   2  53  0.      0.0181  0.      1.025  0.;
   3   4  0.0013  0.0213  0.2214  0.     0.;
   3  18  0.0011  0.0133  0.2138  0.     0.;
   4   5  0.0008  0.0128  0.1342  0.     0.;
   4  14  0.0008  0.0129  0.1382  0.     0.;
   5   6  0.0002  0.0026  0.0434  0.     0.;
   5   8  0.0008  0.0112  0.1476  0.     0.;
   6   7  0.0006  0.0092  0.1130  0.     0.;
   6  11  0.0007  0.0082  0.1389  0.     0.;
   6  54  0.      0.0250  0.      1.07   0.;
   7   8  0.0004  0.0046  0.078   0.     0.;
   8   9  0.0023  0.0363  0.3804  0.     0.;
   9  30  0.0019  0.0183  0.29    0.     0.;
   10  11  0.0004  0.0043  0.0729  0.     0.;
   10  13  0.0004  0.0043  0.0729  0.     0.;
   10  55  0.      0.02    0.      1.07   0.;
   12  11  0.0016  0.0435  0.      1.06   0.;
   12  13  0.0016  0.0435  0.      1.06   0.;
   13  14  0.0009  0.0101  0.1723  0.     0.;
   14  15  0.0018  0.0217  0.366   0.     0.;
   15  16  0.0009  0.0094  0.171   0.     0.;
   16  17  0.0007  0.0089  0.1342  0.     0.;
   16  19  0.0016  0.0195  0.3040  0.     0.;
   16  21  0.0008  0.0135  0.2548  0.     0.;
   16  24  0.0003  0.0059  0.0680  0.     0.;
   17  18  0.0007  0.0082  0.1319  0.     0.;
   17  27  0.0013  0.0173  0.3216  0.     0.;
   19  20  0.0007  0.0138  0.      1.06   0.;
   19  56  0.0007  0.0142  0.      1.07   0.;
   20  57  0.0009  0.0180  0.      1.009  0.;
   21  22  0.0008  0.0140  0.2565  0.     0.;
   22  23  0.0006  0.0096  0.1846  0.     0.;
   22  58  0.      0.0143  0.      1.025  0.;
   23  24  0.0022  0.0350  0.3610  0.     0.;
   23  59  0.0005  0.0272  0.      0.     0.;
   25  26  0.0032  0.0323  0.5310  0.     0.;
   25  60  0.0006  0.0232  0.      1.025  0.;
   26  27  0.0014  0.0147  0.2396  0.     0.;
   26  28  0.0043  0.0474  0.7802  0.     0.;
   26  29  0.0057  0.0625  1.0290  0.     0.;
   28  29  0.0014  0.0151  0.2490  0.     0.;
   29  61  0.0008  0.0156  0.      1.025  0.;
   9  30  0.0019  0.0183  0.29    0.     0.;
   9  36  0.0022  0.0196  0.34    0.     0.;
   9  36  0.0022  0.0196  0.34    0.     0.;
   36  37  0.0005  0.0045  0.32    0.     0.;
   34  36  0.0033  0.0111  1.45    0.     0.;
   35  34  0.0001  0.0074  0.      0.946  0.;
   33  34  0.0011  0.0157  0.202   0.     0.;
   32  33  0.0008  0.0099  0.168   0.     0.;
   30  31  0.0013  0.0187  0.333   0.     0.;
   30  32  0.0024  0.0288  0.488   0.     0.;
   1  31  0.0016  0.0163  0.25    0.     0.;
   31  38  0.0011  0.0147  0.247   0.     0.;
   33  38  0.0036  0.0444  0.693   0.     0.;
   38  46  0.0022  0.0284  0.43    0.     0.;
   46  49  0.0018  0.0274  0.27    0.     0.;
   1  47  0.0013  0.0188  1.31    0.     0.;
   47  48  0.0025  0.0268  0.40    0.     0.;
   47  48  0.0025  0.0268  0.40    0.     0.;
   48  40  0.0020  0.022   1.28    0.     0.;
   35  45  0.0007  0.0175  1.39    0.     0.;
   37  43  0.0005  0.0276  0.      0.     0.;
   43  44  0.0001  0.0011  0.      0.     0.;
   44  45  0.0025  0.073   0.      0.     0.;
   39  44  0.      0.0411  0.      0.     0.;
   39  45  0.      0.0839  0.      0.     0.;
   45  51  0.0004  0.0105  0.72    0.     0.;
   50  52  0.0012  0.0288  2.06    0.     0.;
   50  51  0.0009  0.0221  1.62    0.     0.;
   49  52  0.0076  0.1141  1.16    0.     0.;
   52  42  0.0040  0.0600  2.25    0.     0.;
   42  41  0.0040  0.0600  2.25    0.     0.;
   41  40  0.0060  0.0840  3.15    0.     0.;
   31  62  0.      0.026   0.      1.04   0.;
   32  63  0.      0.013   0.      1.04   0.;
   36  64  0.      0.0075  0.      1.04   0.;
   37  65  0.      0.0033  0.      1.04   0.;
   41  66  0.      0.0015  0.      1.     0.;
   42  67  0.      0.0015  0.      1.     0.;
   52  68  0.      0.0030  0.      1.     0.;
   1  27  0.032   0.32    0.41    1.     0.];



% Machine data format
%       1. machine number,
%       2. bus number,
%       3. base mva,
%       4. leakage reactance x_l(pu),
%       5. resistance r_a(pu),
%       6. d-axis sychronous reactance x_d(pu),
%       7. d-axis transient reactance x'_d(pu),
%       8. d-axis subtransient reactance x"_d(pu),
%       9. d-axis open-circuit time constant T'_do(sec),
%      10. d-axis open-circuit subtransient time constant
%                T"_do(sec),
%      11. q-axis sychronous reactance x_q(pu),
%      12. q-axis transient reactance x'_q(pu),
%      13. q-axis subtransient reactance x"_q(pu),
%      14. q-axis open-circuit time constant T'_qo(sec),
%      15. q-axis open circuit subtransient time constant
%                T"_qo(sec),
%      16. inertia constant H(sec),
%      17. damping coefficient d_o(pu),
%      18. dampling coefficient d_1(pu),
%      19. generator type 
%				1 - classical
%               2 - salient pole
%               3 - round rotor
%      20. saturation factor S(1.0)
%      21. saturation factor S(1.2)
% note: all the following machines use subtransient reactance model

 mac_con = [...
      1  53     300  0.003  	0  	0.969   0.248     0.147  	12.6     0.045  ...
      								0.600   0.250     0      	0.035    0  ...
                                    3.4     0   	  0  		3  	  0.0654  0.5743;% hydro unit
      2  54     800  0.035  	0  	1.8     0.42529   0.30508   6.56  0.05 ...
         							1.7207  0.3661    0.30508   1.5   0.035  ...
                                    4.9494  0   	  0  		3  	  0.0654  0.5743;
      3  55     800	 0.0304  0  	1.8  	0.38309   0.32465   5.7   0.05 ...
      								1.7098  0.36072   0.32465   1.5   0.035  ...
      								4.9623  0   	  0  		3  	  0.0654  0.5743;
   	  4  56  	800  0.0295  0  	1.8  	0.29954   0.24046   5.69  0.05 ...
      								1.7725  0.27481   0.24046   1.5   0.035  ...
      								4.1629  0   	  0  		3  	  0.0654  0.5743;
   	  5  57  	700   0.027  	0  	1.8     0.36  	  0.27273   5.4   0.05 ...
      								1.6909  0.32727   0.27273   0.44  0.035  ...
      								4.7667  0  	      0  		3  	  0.0654  0.5743;
   	  6  58  	900  0.0224  	0  	1.8  	0.35433   0.28346   7.3   0.05 ...
      								1.7079  0.3189    0.28346   0.4   0.035  ...
      								4.9107  0   	  0  		3  	  0.0654  0.5743;
   	 7  59  	800  0.0322  	0  	1.8  	0.29898   0.24407   5.66  0.05 ...
      								1.7817  0.27458   0.24407   1.5   0.035  ...
      								4.3267  0   	  0  		3  	  0.0654  0.5743;
   	 8  60  	800  0.028      0  	1.8  	0.35379   0.27931   6.7   0.05 ...
      								1.7379  0.31034   0.27931   0.41  0.035   ...
      								3.915  	0   	  0  		3  	  0.0654  0.5743;
   	9   61     1000  0.0298  	0  	1.8  	0.48718   0.38462   4.79  0.05 ...
      								1.7521  0.42735   0.38462   1.96  0.035  ...
      								4.0365  0   	  0  		3  	  0.0654  0.5743;
   	10  62     1200  0.0199  	0  	1.8  	0.48675   0.42604   9.37  0.05 ...
      								1.2249  0.47929   0.42604   1.5   0.035  ...
      								2.9106  0   	  0  		3  	  0.0654  0.5743;
   	11  63     1600  0.0103  	0  	1.8  	0.25312   0.16875   4.1   0.05 ...
      								1.7297  0.21094   0.16875   1.5   0.035  ...
      								2.0053  0   	  0  		3  	  0.0654  0.5743;
   	12  64     1900  0.022  	0  	1.8  	0.55248   0.44554   7.4   0.05 ...
      								1.6931  0.49901   0.44554   1.5  	0.035  ...
      								5.1791  0   	  0  		3  	  0.0654  0.5743;
   	13  65    12000  0.003  	0  	1.8  	0.33446   0.24324   5.9   0.05 ...
      								1.7392  0.30405   0.24324   1.5     0.035  ...
      								4.0782  4.0782    0  		  3  	0.0654  0.5743;
   	14  66    10000  0.0017  	0  	1.8    	0.285     0.23       4.1    0.05   ...
      								1.73    0.25      0.23       1.5    0.035       ...
      								3  		3   	  0  		  3  	0.0654  0.5743;
   	15  67    10000  0.0017  	0  	1.8    	0.285     0.23       4.1    0.05   ...
      								1.73    0.25      0.23   	 1.5    0.035       ...
      								3  		3   	  0  		 3  	0.0654  0.5743;
   	16  68    11000  0.0041  	0  	1.8  	0.35899   0.27809    7.8    0.05 ...
      								1.6888  0.30337   0.27809    1.5    0.035    ...
                                    4.45  	4.45      0  		 3  	0.0654  0.5743;
 ];
% exciter data dc12 model
%     1 - exciter type (1 for DC1, 2 for DC2)
%     2 - machine number
%     3 - input filter time constant T_R
%     4 - voltage regulator gain K_A
%     5 - voltage regulator time constant T_A
%     6 - voltage regulator time constant T_B
%     7 - voltage regulator time constant T_C
%     8 - maximum voltage regulator output V_Rmax
%     9 - minimum voltage regulator output V_Rmin
%    10 - exciter constant K_E
%    11 - exciter time constant T_E
%    12 - E_1
%    13 - saturation function S_E(E_1)
%    14 - E_2
%    15 - saturation function S_E(E_2)
%    16 - stabilizer gain K_F
%    17 - stabilizer time constant T_F

exc_condc = [...
   1 1 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 2 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 3 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 4 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 5 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 6 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 7 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 8 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 9 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0
   1 10 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 11 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0;
   1 12 0 30. 0.02 0 0 10.  -10.  1.0    .785 3.9267 0.070 5.2356 0.910 0.030 1.0 0 0 0];

exc_con = [...
0 1 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 2 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 3 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 4 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 5 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 6 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 7 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 8 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 9 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 10 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 11 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 12 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0; 
0 13 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 14 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 15 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0; 
0 16 0.0 100.0  0.01     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;  
    ];
pss_con = [...
     1   1    100  10  0.1   0.02  0.08  0.02 0.2 -0.05 ; 
     1   2    100  10  0.08  0.02  0.08  0.02 0.2 -0.05 ;   
     1   3    100  10  0.08  0.02  0.08  0.02 0.2 -0.05 ;
     1   4    100  10  0.08  0.02  0.08  0.02 0.2 -0.05 ;
     1   5    100  10  0.08  0.02  0.08  0.02 0.2 -0.05 ;
     1   6    100  10  0.1   0.02  0.1   0.02 0.2 -0.05 ;
     1   7    100  10  0.08  0.02  0.08  0.02 0.2 -0.05 ;
     1   8    100  10  0.08  0.02  0.08  0.02 0.2 -0.05 ;
     1   9    100  10  0.08  0.03  0.05  0.01 0.2 -0.05 ; 
     1   10   100  10  0.1   0.02  0.1   0.02 0.2 -0.05 ;
     1   11    50  10  0.08  0.03 0.05   0.01 0.2 -0.05 ;
     1   12   110  10  0.1   0.02  0.1   0.02 0.2 -0.05 ;
];
%pss_con = [];
% non-conforming load on all load buses
load_con = [...
   1  0 0 .5 0;
   3  0 0 .5 0;
   4  0 0 .5 0;
   7  0 0 .5 0;
   8  0 0 .5 0;
   9  0 0 .5 0;
   15 0 0 .5 0;
   16 0 0 .5 0;
   18 0 0 .5 0;
   20 0 0 .5 0;
   21 0 0 .5 0;
   23 0 0 .5 0;
   24 0 0 .5 0;
   25 0 0 .5 0;
   26 0 0 .5 0;
   27 0 0 .5 0;
   28 0 0 .5 0;
   29 0 0 .5 0;
   33 0 0 .5 0;
   37 0 0 .5 0;
   39 0 0 .5 0;
   40 0 0 .5 0;
   41 0 0 .5 0;
   42 0 0 .5 0;
   44 0 0 .5 0;
   45 0 0 .5 0;
   46 0 0 .5 0;
   47 0 0 .5 0;
   48 0 0 .5 0;
   49 0 0 .5 0;
   50 0 0 .5 0;
   51 0 0 .5 0;
   52 0 0 .5 0];
ind_con = [ ...
  1   1  150  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  2   3  200  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  3   4  350  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  4   7  150  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  5   8  325  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  6   9   60  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  7   12   5  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  8   15 200  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  9   16 200  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  10  18 100  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  11  20 400  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  12  21 100  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  13  23 100  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  14  24 175  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  15  25 150  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  16  26 100  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  17  27 200  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  18  28 120  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  19  29 200  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  20  33  75  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  21  36  75  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  22  40  50  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  23  45 120  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  24  46 100  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  25  47 120  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  26  48 200  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  27  49 110  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  28  50 130  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
  29  51 300  .03274 .08516 3.7788  .06164 .06005  1.0 0.01354 0.07517 0     3 0 .5;
];
mld_con = [ ...
1  1  .1  1  .7  2;
2  3  .1  1  .7  2;
3  4  .1  1  .7  2;
4  7  .1  1  .7  2;
5  8  .1  1  .7  2;
6  9  .1  1  .7  2;
7  12 .1  1  .7  2;
8  15 .1  1  .7  2;
9  16 .1  1  .7  2;
10 18 .1  1  .7  2;
11 20 .1  1  .7  2;
12 21 .1  1  .7  2;
13 23 .1  1  .7  2;
14 24 .1  1  .7  2;
15 25 .1  1  .7  2;
16 26 .1  1  .7  2;
17 27 .1  1  .7  2;
18 28 .1  1  .7  2;
19 29 .1  1  .7  2;
20 33 .1  1  .7  2;
21 36 .1  1  .7  2;
22 40 .1  1  .7  2;
23 45 .1  1  .7  2;
24 46 .1  1  .7  2;
25 47 .1  1  .7  2;
26 48 .1  1  .7  2;
27 49 .1  1  .7  2;
28 50 .1  1  .7  2;
29 51 .1  1  .7  2;
];

%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault  - 0 three phase
%                            - 1 line to ground
%                            - 2 line-to-line to ground
%                            - 3 line-to-line
%                            - 4 loss of line with no fault
%                            - 5 loss of load at bus
%       col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)
%
%
%
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)

sw_con = [...
   0    0    0    0    0    0    0.01;%sets intitial time step
   0.1  1    2   0    0    0    0.005; %apply three phase fault at bus 1, on line 1-2
   0.15 0    0    0    0    0    0.005; %clear fault at bus 
   0.20 0    0    0    0    0    0.005; %clear remote end
   0.50 0    0    0    0    0    0.01; % increase time step 
   1.0  0    0    0    0    0    0.01; % increase time step
   10.0  0    0    0    0    0    0]; % end simulation
%
%
% non-conforming load
% col 1           bus number
% col 2           fraction const active power load
% col 3           fraction const reactive power load
% col 4           fraction const active current load
% col 5           fraction const reactive current load

