function out = interpLsp(p,Rc)

% Does interpolation in Marrero and others to get subsurface effective
% attenuation length for given atmospheric pressure (hPa) and cutoff
% rigidity (GV). 
%
% out = interpLsp(p,Rc)
% 
% Vectorized. Send vectors of common length.


pp = [400         500         600         700         800         900        1000        1100
         400         500         600         700         800         900        1000        1100
         400         500         600         700         800         900        1000        1100
         400         500         600         700         800         900        1000        1100
         400         500         600         700         800         900        1000        1100
         400         500         600         700         800         900        1000        1100];
     
RcRc = [0     0     0     0     0     0     0     0
     4     4     4     4     4     4     4     4
     8     8     8     8     8     8     8     8
    12    12    12    12    12    12    12    12
    16    16    16    16    16    16    16    16
    20    20    20    20    20    20    20    20];

LL  = [154   146   142   139   138   137   136   136
   166   154   147   143   140   138   137   136
   183   166   157   151   147   144   142   139
   195   176   165   158   154   150   148   145
   205   183   171   164   159   155   152   151
   207   185   172   165   160   156   153   151];

if any(p < 400);
    p(find(p < 400)) = 400 + zeros(size(find(p < 400)));
end;

if any(p > 1100);
    p(find(p > 1100)) = 1100 + zeros(size(find(p > 1100)));
end;

if any(Rc > 20);
    Rc(find(Rc > 20)) = 20 + zeros(size(find(Rc > 20)));
end;


out = interp2(pp,RcRc,LL,p,Rc,'linear');
