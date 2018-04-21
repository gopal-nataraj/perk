function [X_flip] = flipud3D(X)
%flipud3D Slice-by-slice flipud() for 3D inputs

[nx, ny, nz] = size(X);
X_flip = NaN(nx, ny, nz);
for z = 1:nz
    X_flip(:,:,z) = flipud(squeeze(X(:,:,z)));
end

end

