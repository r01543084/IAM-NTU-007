function [ux0,uy0,T0,density0,p0] = IC(configuration,nx,ny)
switch configuration
    case(1)
        p       = [1        0.4         0.0439      0.15    ];%pressure
        density = [1        0.5197      0.1072      0.2579  ];%density
        ux      = [0        -0.7259     -0.7259     0       ];%x velocity
        uy      = [0        0           -1.4045     -1.4045	];%y velocity
        T = p./(density);    % Temperature
	case(2)
        p       = [1        0.4         1           0.4     ];%pressure
        density = [1        0.5197      1           0.5197  ];%density
        ux      = [0        -0.7259     -0.7259     0       ];%x velocity
        uy      = [0        0           -0.7259     -0.7259	];%y velocity
        T = p./(density);    % Temperature
	case(3)
        p       = [1.5      0.3         0.029       0.3     ];%pressure
        density = [1.5      0.5232      0.138       0.5323  ];%density
        ux      = [0        1.206       1.206       0       ];%x velocity
        uy      = [0        0           1.206       1.206	];%y velocity
        T = p./(density);    % Temperature
	case(4)
        p       = [1.1      0.35        1.1         0.35    ];%pressure
        density = [1.1      0.5065      1.1         0.5065  ];%density
        ux      = [0        0.8939      0.8939      0       ];%x velocity
        uy      = [0        0           0.8939      0.8939	];%y velocity
        T = p./(density);    % Temperature
	case(5)
        p       = [1        1           1           1       ];%pressure
        density = [1        2           1           3       ];%density
        ux      = [-0.75    -0.75       0.75        0.75    ];%x velocity
        uy      = [-0.5     0.5         0.5         -0.5    ];%y velocity
        T = p./(density);    % Temperature
	case(6)
        p       = [1        1           1           1       ];%pressure
        density = [1        2           1           3       ];%density
        ux      = [0.75     0.75        -0.75       -0.75   ];%x velocity
        uy      = [-0.5     0.5         0.5         -0.5    ];%y velocity
        T = p./(density);    % Temperature	
    case(7)
        p       = [1        0.4         0.4         0.4     ];%pressure
        density = [1        0.5197      0.8         0.5197	];%density
        ux      = [0.1      -0.6259     0.1         0.1     ];%x velocity
        uy      = [0.1      0.1         0.1         -0.6259	];%y velocity
        T = p./(density);    % Temperature
    case(8)
        p       = [0.4      1           1           1       ];%pressure
        density = [0.5197   1           0.8         1   	];%density
        ux      = [0.1      -0.6259     0.1         0.1     ];%x velocity
        uy      = [0.1      0.1         0.1         -0.6259	];%y velocity
        T = p./(density);    % Temperature
	case(9)
        p       = [1        1           0.4         0.4     ];%pressure
        density = [1        2           1.039       0.5197  ];%density
        ux      = [0        0           0           0       ];%x velocity
        uy      = [0.3      -0.3        -0.8133     -0.4259	];%y velocity
        T = p./(density);    % Temperature
	case(10)
        p       = [1        1           0.3333		0.3333	];%pressure
        density = [1        0.5     	0.2281     	0.4562  ];%density
        ux      = [0        0           0           0       ];%x velocity
        uy      = [0.4297	0.6076    	-0.6076     -0.4297	];%y velocity
        T = p./(density);    % Temperature
	case(11)
        p       = [1        0.4         0.4 		0.4 	];%pressure
        density = [1        0.5313     	0.8      	0.5313  ];%density
        ux      = [0.1   	0.8276      0.1         0.1     ];%x velocity
        uy      = [0        0           0           0.7276	];%y velocity
        T = p./(density);    % Temperature
	case(12)
        p       = [0.4      1           1           1       ];%pressure
        density = [0.5313   1           0.8      	1       ];%density
        ux      = [0        0.7276      0           0       ];%x velocity
        uy      = [0        0           0           0.7276	];%y velocity
        T = p./(density);    % Temperature
	case(13)
        p       = [1        1           0.4         0.4     ];%pressure
        density = [1        2           1.0625      0.5313  ];%density
        ux      = [0        0           0           0       ];%x velocity
        uy      = [-0.3     0.3         0.8145      0.4276	];%y velocity
        T = p./(density);    % Temperature
	case(14)
        p       = [8        8           2.6667      2.6667  ];%pressure
        density = [2        1           0.4736      0.9474  ];%density
        ux      = [0        0           0           0       ];%x velocity
        uy      = [-0.5606	-1.2172     1.2172      1.1606	];%y velocity
        T = p./(density);    % Temperature
	case(15)
        p       = [1        0.4         0.4         0.4    ];%pressure
        density = [1        0.5197      0.8         0.5313  ];%density
        ux      = [0.1      -0.6259     0.1         0.1     ];%x velocity
        uy      = [-0.3     -0.3        -0.3        0.4276	];%y velocity
        T = p./(density);    % Temperature
	case(16)
        p       = [0.4      1           1           1       ];%pressure
        density = [0.5313   1.0222      0.8         1       ];%density
        ux      = [0.1      -0.6179     0.1         0.1     ];%x velocity
        uy      = [0.1      0.1         0.1         0.8276	];%y velocity
        T = p./(density);    % Temperature
	case(17)
        p       = [1        1           0.4         0.4     ];%pressure
        density = [1        2           1.0625      0.5197  ];%density
        ux      = [0        0           0           0       ];%x velocity
        uy      = [-0.4     -0.3        0.2145      -1.1259 ];%y velocity
        T = p./(density);    % Temperature
	case(18)
        p       = [1        1           0.4         0.4     ];%pressure
        density = [1        2           1.0625      0.5197  ];%density
        ux      = [0        0           0           0       ];%x velocity
        uy      = [1        -0.3        0.2145      0.2741  ];%y velocity
        T = p./(density);    % Temperature
	case(19)
        p       = [1        1           0.4         0.4     ];%pressure
        density = [1        2           1.0625      0.5197  ];%density
        ux      = [0        0           0           0       ];%x velocity
        uy      = [0.3      -0.3        0.2145      -0.4259 ];%y velocity
        T = p./(density);    % Temperature
	case(20)
        p       = [0.1      0.5         0.5         0.1     ];%pressure
        density = [0.1      7.3/3       7.3/3       0.1     ];%density
        ux      = [0        1.1711      1.1711      0       ];%x velocity
        uy      = [0        0           0           0       ];%y velocity
        T = p./(density);    % Temperature
end


    
    %¨ú¤¤ÂI
%     x_middle = ceil(nx/2);
%     y_middle = ceil(ny/2);
%     
%     stagex1 = 1:x_middle; stagex2 = x_middle+1:nx;
%     stagey1 = 1:y_middle; stagey2 = y_middle+1:ny;
    stagex1 = 1:10; stagex2 = 11:nx;
    stagey1 = 1:10; stagey2 = 11:ny;
    
    % Initial Condition for our 2D domain
    % Velovity in x
    ux0(stagey2,stagex2) = ux(1); % region 1
    ux0(stagey2,stagex1) = ux(2); % region 2
    ux0(stagey1,stagex1) = ux(3); % region 3
    ux0(stagey1,stagex2) = ux(4); % region 4
	% Velovity in y
    uy0(stagey2,stagex2) = uy(1); % region 1
    uy0(stagey2,stagex1) = uy(2); % region 2
    uy0(stagey1,stagex1) = uy(3); % region 3
    uy0(stagey1,stagex2) = uy(4); % region 4
    % Temperature
    T0(stagey2,stagex2) = T(1); % region 1
    T0(stagey2,stagex1) = T(2); % region 2
    T0(stagey1,stagex1) = T(3); % region 3
    T0(stagey1,stagex2) = T(4); % region 4
    % density
    density0(stagey2,stagex2) = density(1); % region 1
    density0(stagey2,stagex1) = density(2); % region 2
    density0(stagey1,stagex1) = density(3); % region 3
    density0(stagey1,stagex2) = density(4); % region 4
    % pressure
    p0(stagey2,stagex2) = p(1); % region 1
    p0(stagey2,stagex1) = p(2); % region 2
    p0(stagey1,stagex1) = p(3); % region 3
    p0(stagey1,stagex2) = p(4); % region 4
    