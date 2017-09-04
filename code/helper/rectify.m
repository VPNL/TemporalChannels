function X_rect = rectify(X)
% Rectify response time series X by setting negative values to zero.
% AS 9/2017

X_rect = X;
X_rect(X < 0) = 0;

end
