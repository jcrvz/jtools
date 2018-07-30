% function A = pivot(A,r,s)
%
% ars = A(r,s);
%
% % Row operations
% for ir = 1 : size(A,1)
%     if ir ~= r
%
%         airs = A(ir,s);
%         A(ir,:) = ars*A(ir,:) - airs*A(r,:);
%     end
% end
%
% % Making ones be ones
% for is = 1 : size(A,2)
%     elements = find(A(:,is) ~= 0);
%     if numel(elements) == 1
%         A(elements,:) = A(elements,:)/A(elements,is);
%     end
% end

function R = pivot(M, r, c)
[d, w] = size(M);             % Get matrix dimensions
R = zeros(d, w);              % Initialize to appropriate size
R(r,:) = M(r, :) / M(r,c);    % Copy row r, normalizing M(r,c) to 1
for k = 1:d                   % For all matrix rows
    if (k ~= r)                 % Other then r
        R(k,:) = M(k,:) ...       % Set them equal to the original matrix
            - M(k,c) * R(r,:); % Minus a multiple of normalized row r, making R(k,c)=0
    end
end
end