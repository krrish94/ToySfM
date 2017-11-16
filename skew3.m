function omega_cross = skew3(omega)
% SKEW3  Takes in a 3-vector omega and constructs a 3 x 3 skew-symmetric
% matrix omega_cross, such that omega_cross * b is same as cross(omega, b).

omega_cross = [0, -omega(3), omega(2); ...
                omega(3), 0, -omega(1); ...
                -omega(2), omega(1), 0];

end
