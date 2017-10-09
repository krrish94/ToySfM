function pts = generateCube(n, a)
% GENERATECUBE  Synthesize a cube (3-D) containing n-by-n-by-n points and
% having side length 'a'.

% Determine the distance between successive points on the cube
% Note the use of 'n-1' as opposed to 'n'. If we divide a into n-1 equal
% intervals, we get n points.
delta = a / (n-1);

% Generate points on the cube
% pts = zeros(3, n*n*n);
pts = [];
idx = 1;
for i = -a/2:delta:a/2
    for j = -a/2:delta:a/2
        for k = -a/2:delta:a/2
            % pts(:,idx) = [i; j; k];
            % Accept a point only if it is on the surface of the cube
            if abs(i) == a/2 || abs(j) == a/2 || abs(k) == a/2
            % if max([abs(i), abs(j), abs(k)]) == a/2
                pts = [pts, [i; j; k]];
            end
            % idx = idx + 1;
        end
    end
end

end
