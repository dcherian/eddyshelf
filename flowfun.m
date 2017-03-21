function  [phi,psi] = flowfun(xv, yv, u,v,flag)

% FLOWFUN  Computes the potential PHI and the streamfunction PSI
%     of a 2-dimensional flow defined by the matrices of velocity
%     components U and V, so that
%
%           d(PHI)    d(PSI)          d(PHI)    d(PSI)
%      u =  -----  -  ----- ,    v =  -----  +  -----
%            dx        dy              dy        dx
%
%     For a potential (irrotational) flow  PSI = 0, and the laplacian
%     of PSI is equal to the divergence of the velocity field.
%     A non-divergent flow can be described by the streamfunction
%     alone, and the laplacian of the streamfunction is equal to
%     vorticity (curl) of the velocity field.
%     The stepsizes dx and dy are assumed to equal unity for PHI.
%   [PHI,PSI] = FLOWFUN(U,V), or in a complex form
%   [PHI,PSI] = FLOWFUN(U+iV)
%     returns matrices PHI and PSI of the same sizes as U and V,
%     containing potential and streamfunction given by velocity
%     components U, V.
%     Because these potentials are defined up to the integration
%     constant their absolute values are such that
%     PHI(1,1) = PSI(1,1) = 0.
%     If only streamfunction is needed, the flag can be used:
%   PSI = FLOWFUN(U,V,FLAG), where FLAG can be a string:
%     '-', 'psi', 'streamfunction' (abbreviations allowed).
%     For the potential the FLAG can be  '+', 'phi', 'potential'.

%  Uses command CUMSIMP (Simpson rule summation).

%  Kirill K. Pankratov, March 7, 1994.

% Check input arguments .............................................
issu=0; issv=0; isflag=0;    % For input checking
isphi = 1; ispsi = 1;        % Is phi and psi to be computed
if nargin==1, issu = ischar(u); end
if nargin==2, issv = ischar(v); end
if nargin==1 && ~issu, v=imag(u); end
if issv, flag = v; v = imag(u); isflag = 1; end 
if nargin==0 || issu            % Not enough input arguments
  disp([10 '  Error: function must have input arguments:'...
  10 '  matrivces  U and V  (or complex matrix W = U+iV)' 10 ])
  return
end
if any(size(u)~=size(v))     % Disparate sizes
  disp([10 '  Error: matrices U and V must be of equal size' 10])
  return
end
if nargin==3, isflag=1; end
u = real(u);

 % Check the flag string . . . . . . . .
Dcn = char('+','potential','phi');
Dcn = char(Dcn,'-','streamfunction','psi');
if isflag
  lmin = min(size(flag,2),size(Dcn,2));
  flag = flag(1,1:lmin);
  A = flag(ones(size(Dcn,1),1),1:lmin)==Dcn(:,1:lmin);
  if lmin>1, coinc = sum(A'); else coinc = A'; end
  fnd = find(coinc==lmin);
  if ~isempty(fnd), if fnd<4, ispsi=0; else isphi=0; end, end
end

% Now the main computations .........................................
% Integrate velocity fields to get potential and streamfunction
% Use Simpson rule summation (function CUMSIMP)

 % Compute potential PHI (potential, non-rotating part)
if isphi
  cx = cumtrapz(u(:,1,:),1);  % Compute x-integration constant
  cy = cumtrapz(v(1,:,:),2);  % Compute y-integration constant
  phi = bsxfun(@plus, cumtrapz(yv, v, 2), cx);
  phi = bsxfun(@plus, phi+cumtrapz(xv, u, 1), cy)/2;
end

 % Compute streamfunction PSI (solenoidal part)
if ispsi
  cx = cumtrapz(xv, v(:,1,:),1);  % Compute x-integration constant
  cy = cumtrapz(yv, u(1,:,:),2);  % Compute y-integration constant
  psi = bsxfun(@plus, -cumtrapz(yv, u,2), cx);
  psi = bsxfun(@plus, psi+cumtrapz(xv, v,1), -cy)/2;
end

 % Rename output if need only PSI
if ~isphi && ispsi && nargout==1, phi = psi; end

end
%=========================== end  flowfun.m ======================

%============================ save as  cumsimp.m =================
function  f = cumsimp(y)

% F = CUMSIMP(Y)    Simpson-rule column-wise cumulative summation.
%       Numerical approximation of a function F(x) such that 
%       Y(X) = dF/dX.  Each column of the input matrix Y represents
%       the value of the integrand  Y(X)  at equally spaced points
%       X = 0,1,...size(Y,1).
%       The output is a matrix  F of the same size as Y.
%       The first row of F is equal to zero and each following row
%       is the approximation of the integral of each column of matrix
%       Y up to the givem row.
%       CUMSIMP assumes continuity of each column of the function Y(X)
%       and uses Simpson rule summation.
%       Similar to the command F = CUMSUM(Y), exept for zero first
%       row and more accurate summation (under the assumption of
%       continuous integrand Y(X)).
% 
%    See also CUMSUM, SUM, TRAPZ, QUAD

%  Kirill K. Pankratov, March 7, 1994.

 % 3-points interpolation coefficients to midpoints.
 % Second-order polynomial (parabolic) interpolation coefficients
 % from  Xbasis = [0 1 2]  to  Xint = [.5 1.5]
c1 = 3/8; c2 = 6/8; c3 = -1/8;

 % Determine the size of the input and make column if vector
ist = 0;         % If to be transposed
lv = size(y,1);
if lv==1, ist = 1; y = y(:); lv = length(y); end
f = zeros(size(y));

 % If only 2 elements in columns - simple sum divided by 2
if lv==2
  f(2,:) = (y(1,:)+y(2))/2;
  if ist, f = f'; end   % Transpose output if necessary
  return
end

 % If more than two elements in columns - Simpson summation
num = 1:lv-2;
   % Interpolate values of Y to all midpoints
f(num+1,:) = c1*y(num,:)+c2*y(num+1,:)+c3*y(num+2,:);
f(num+2,:) = f(num+2,:)+c3*y(num,:)+c2*y(num+1,:)+c1*y(num+2,:);
f(2,:) = f(2,:)*2; f(lv,:) = f(lv,:)*2;
   % Now Simpson (1,4,1) rule
f(2:lv,:) = 2*f(2:lv,:)+y(1:lv-1,:)+y(2:lv,:);
f = cumsum(f)/6;  % Cumulative sum, 6 - denom. from the Simpson rule

if ist, f = f'; end     % Transpose output if necessary

end
%============================= end  cumsimp.m =================
