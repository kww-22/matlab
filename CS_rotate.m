GtL = [.926 -.052 .374; -.013 .985 .170; -.377 -.162 .912];

LtG = GtL';

% actual = GtL*roty(-90);
% 
% correct = actual*LtG;
% 
% disp(correct)

sensor = [.995 -.100 -.000; .099 .981 .169; -.017 -.168 .986];

correct = sensor*LtG;

disp(correct)