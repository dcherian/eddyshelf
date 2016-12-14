
vsamp = 4/2; % m/s
vsteam = 4/2; % m/s
vsteamas = vsamp;

80e3/vsamp/3600 + ... %shelfbreak along-shelf line
    sqrt(2)*20e3/vsamp/3600 + ... % diagonal to offshore line 1
    70e3/vsamp/3600 + ... % offshore along-shelf line 1
    20e3/vsamp/3600 + ... % diagonal to offshore line 1
    70e3/vsamp/3600 + ... % offshore along-shelf line 2
    sqrt(2)*20e3/vsamp/3600 + ... % diagonal to cross-shelf line 1
    60e3/vsamp/3600 + ... % cross-shelf line 1
    20e3/vsteamas/3600 + ... % along-shelf steam
    60e3/vsamp/3600 + ... % cross-shelf line 2
    20e3/vsteamas/3600 + ... % along-shelf steam
    60e3/vsamp/3600 + ... % cross-shelf line 3
    20e3/vsteamas/3600 + ... % along-shelf steam
    60e3/vsamp/3600 + ... % cross-shelf line 4
    20e3/vsteamas/3600 + ... % along-shelf steam
    60e3/vsamp/3600 + ... % cross-shelf line 5
    40e3/vsamp/3600 % return on 5 to shelfbreak