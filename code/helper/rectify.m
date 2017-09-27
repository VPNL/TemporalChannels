function X_rect = rectify(X, polarity)
% Rectify data (X) in the positive or negative direction (polarity).
% AS 9/2017

if nargin < 2
    polarity = 'positive';
end
X_rect = X;
switch polarity
    case 'positive'
        X_rect(X < 0) = 0;
    case 'negative'
        X_rect(X > 0) = 0;
    otherwise
        error('Polarity must be positive or negative');
end

end
