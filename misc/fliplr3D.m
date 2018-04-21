function [X_flip] = fliplr3D(X)
%fliplr3D Slice-by-slice fliplr() for 3D inputs

[nx, ny, nz] = size(X);
X_flip = NaN(nx, ny, nz);
for z = 1:nz
    X_flip(:,:,z) = fliplr(squeeze(X(:,:,z)));
end

end

