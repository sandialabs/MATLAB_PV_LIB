classdef GHIImage
    properties
        Value
    end
    
    methods
        function obj = GHIImage(val)
            if nargin>0
                if ismatrix(val) && ~isvector(val) &&(isnumeric(val) | islogical(val))
                    obj.Value = val;
                else
                    error('GHIIMage must be 2D numeric or logical array');
                end
            end
        end
                    
        function matrixOut = smooth2a(matrixIn,Nr,Nc)
            % Smooths 2D array data.  Ignores NaN's.
            %
            % function matrixOut = smooth2a(matrixIn,Nr,Nc)
            % 
            % This function smooths the data in matrixIn using a mean filter over a
            % rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
            % element "i" by the mean of the rectange centered on "i".  Any NaN
            % elements are ignored in the averaging.  If element "i" is a NaN, then it
            % will be preserved as NaN in the output.  At the edges of the matrix,
            % where you cannot build a full rectangle, as much of the rectangle that
            % fits on your matrix is used (similar to the default on Matlab's builtin
            % function "smooth").
            % 
            % "matrixIn": original matrix
            % "Nr": number of points used to smooth rows
            % "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
            % 
            % "matrixOut": smoothed version of original matrix
            % 
            % 
            % 	Written by Greg Reeves, March 2009.
            % 	Division of Biology
            % 	Caltech
            % 
            % 	Inspired by "smooth2", written by Kelly Hilands, October 2004
            % 	Applied Research Laboratory
            % 	Penn State University
            % 
            % 	Developed from code written by Olof Liungman, 1997
            % 	Dept. of Oceanography, Earth Sciences Centre
            % 	G�teborg University, Sweden
            % 	E-mail: olof.liungman@oce.gu.se

            %
            % Initial error statements and definitions
            %
            if nargin < 2, error('Not enough input arguments!'), end

            N(1) = Nr; 
            if nargin < 3, N(2) = N(1); else N(2) = Nc; end

            if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
            if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

            %
            % Building matrices that will compute running sums.  The left-matrix, eL,
            % smooths along the rows.  The right-matrix, eR, smooths along the
            % columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by- 
            % (2*Nc+1) rectangle centered on element "i".
            %
            [row,col] = size(matrixIn.Value);
            eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
            eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

            %
            % Setting all "NaN" elements of "matrixIn" to zero so that these will not
            % affect the summation.  (If this isn't done, any sum that includes a NaN
            % will also become NaN.)
            %
            A = isnan(matrixIn.Value);
            matrixIn.Value(A) = 0;

            %
            % For each element, we have to count how many non-NaN elements went into
            % the sums.  This is so we can divide by that number to get a mean.  We use
            % the same matrices to do this (ie, "eL" and "eR").
            %
            nrmlize = eL*(~A)*eR;
            nrmlize(A) = NaN;

            %
            % Actually taking the mean.
            %
            matrixOut = eL*matrixIn.Value*eR;
            matrixOut = matrixOut./nrmlize;
        end
        
        function out = bwareaopen(X,connectedness,threshold)
            % BWAREAOPEN: perform area opening of binary image X
            %
            %   Out = BWAREAOPEN(X, Y,SE)
            %
            %   Inputs:
            %     X: input image array (assumed binary [0 1])
            %
            %   Outputs
            %     Out: output image array
            if ~islogical(X.Value)
                error('bwareaopen passed non-binary array')
            end

            cc = findCC(X,connectedness);  % Get connected components
            counts = histcounts(cc,'BinMethod','integers');
            % integers = 0:65000;
            % counts = hist(cc(:),integers);
            out = ismember(cc,find(counts(2:end)>threshold));
        end
        
        function out = bwclose(img, SE)
            % BWCLOSE: morphological close of binary image with structuring element SE
            %
            %   Out = BWCLOSE(In,SE)
            %
            %   Inputs:
            %     In: input image array (assumed binary [0 1])
            %     SE: structuring element (assumed binary)
            %
            %   Outputs
            %     Out: output image array
            if ~islogical(img.Value)
                error('bwclose passed non-binary array')
            end

            % Version 1: use convolution
            tmp = img.bwdilate(SE);
            out = GHIImage(tmp).bwerode(SE);
        end
        
        function out = bwdilate(img, SE)
            % BWDILATE: Dilate a binary image with structuring element SE
            %
            %   Out = BWDILATE(In,SE)
            %
            %   Inputs:
            %     In: input image array (assumed binary [0 1])
            %     SE: structuring element (assumed binary)
            %
            %   Outputs
            %     Out: output image array
            if ~islogical(img.Value)
                error('bwdilate passed non-binary array')
            end

            % Version 1: use convolution
            out = conv2(single(img.Value), single(SE),'same') > 0;
        end

        function out = bwerode(img, SE)
            % BWERODE: Erode a binary image with structuring element SE
            %
            %   Out = BWERODE(In,SE)
            %
            %   Inputs:
            %     In: input image array (assumed binary [0 1])
            %     SE: structuring element (assumed binary)
            %
            %   Outputs
            %     Out: output image array
            if ~islogical(img.Value)
                error('bwerode passed non-binary array')
            end
            
            % Version 2: erosion of foreground is dilation of background
            out = ~(conv2(single(~img.Value), single(SE), 'same') > 0);

            % Version 1: use convolution
            % out = conv2(single(img),single(SE),'same') == conv2(ones(size(img)),SE,'same');
            % Previously, I used divide, which probably takes longer than array ==
            % out = conv2(single(img),single(SE),'same')./conv2(ones(size(img)),SE,'same') == 1;
 
        end
        
        function out = bwhitmiss(X,J,K)
            % BWHITMISS: perform hit or miss operation on binary image X 
            % 
            %   BWHITMISS(X,J) uses the same structuring element in both phases of
            %   the hit-miss operation: BWHITMISS(X,J) = BWERODE(X,J) & BWERODE(~X,~J) 
            %
            %   BWHITMISS(X,J,K) accommodates "don't care" by allowing different
            %   structuring elements: BWHITMISS(X,J,K) = BWERODE(X,J) & BWERODE(~X,K) 
            %   In this form, K should be the complement of J, with "don't cares" reset
            %   to 0.
            %
            %   Usage: Out = BWHITMISS(X,J) or Out = BWHITMISS(X,J,K)
            %
            %   Inputs:
            %     X: input image array (assumed logical)
            %     J, K: binary structuring elements
            %
            %   Outputs
            %     Out: output (binary/logical) image array

            if ~islogical(X.Value)
                error('bwhitmiss passed non-binary array')
            end

            if nargin < 3
                K = ~J;
            end
            notX = GHIImage(~X.Value);
            out = X.bwerode(J) & notX.bwerode(K);
        end

        function out = bwlengthopen(img,connectedness,threshold)
            % BWAREAOPEN: eliminate connected components in img shorter than threshold
            %
            %   Out = BWAREAOPEN(img,Connectedness,MinLength)
            %
            %   Inputs:
            %     img: input image array (assumed binary [0 1])
            %     Connectedness: 4 or 8
            %     MinLength: minimum length of connected component to retain
            %
            %   Outputs
            %     Out: output image array

            if ~islogical(img.Value)
                error('bwlengthopen passed non-binary array')
            end
            
            if ~ismember([4 8], connectedness)
                error('bwlengthopen: connectedness must be 4 or 8')
            end
            
            cc = findCC(img,connectedness);  % Get connected components

            for ii = 1:max(cc(:))
                [~,cols] = ind2sub(size(cc),find(cc == ii));
                len(ii) = max(cols)-min(cols)+1;
            end

            out = ismember(cc,find(len>threshold));
        end

        function out = bwopen(img, SE)
            % BWOPEN: morphological open of binary image with structuring element SE
            %
            %   Out = BWOPEN(In,SE)
            %
            %   Inputs:
            %     In: input image array (assumed binary [0 1])
            %     SE: structuring element (assumed binary)
            %
            %   Outputs
            %     Out: output image array

            % Version 1: use convolution
            tmp = img.bwerode(SE);
            
            out = GHIImage(tmp).bwdilate(SE);
        end

        function out = graydilate(img, SE)

            % img: grayscale image 
            % SE:  structuring element (assumed binary)

            [M,N] = size(img.Value);
            [SEM, SEN] = size(SE);

            halfSEM = floor(SEM/2);
            halfSEN = floor(SEN/2);

            bufM = M+2*halfSEM;
            bufN = N+2*halfSEN;

            % Create vector of offsets
            offsets = bsxfun(@plus, bufM*(-halfSEN:halfSEN), (-halfSEM:halfSEM)');
            offsets = offsets(SE>0);

            % Put input image into a buffer (-inf guarantees nonselection by max())
            buf = -inf*ones(bufM, bufN);
            buf(halfSEM+(1:M),halfSEN+(1:N)) = img.Value;

            % Iterate through original image pixels
            out = zeros(M, N);
            index = offsets + halfSEN * bufM + halfSEM;
            for jj = 1:N
                for ii = 1:M
                    index = index + 1;
                    out(ii,jj) = max(buf(index));
                end
                % Skip to beginning of next column of original image
                index = index + 2*halfSEM;
            end
        end
        
        function out = grayerode(img, SE)

            % img: grayscale image 
            % SE:  structuring element (assumed binary)

            [M,N] = size(img.Value);
            [SEM, SEN] = size(SE);

            halfSEM = floor(SEM/2);
            halfSEN = floor(SEN/2);

            bufM = M+2*halfSEM;
            bufN = N+2*halfSEN;

            % Create vector of offsets
            offsets = bsxfun(@plus, bufM*(-halfSEN:halfSEN), (-halfSEM:halfSEM)');
            offsets = offsets(SE>0);

            % Put input image into a buffer (inf guarantees nonselection by min())
            buf = inf*ones(bufM, bufN);
            buf(halfSEM+(1:M),halfSEN+(1:N)) = img.Value;

            % Iterate through original image pixels by columns (reduces multiplies inside loop)
            out = zeros(M, N);
            index = offsets + halfSEN * bufM + halfSEM;
            for jj = 1:N
                for ii = 1:M
                    index = index + 1;
                    out(ii,jj) = min(buf(index));
                end
                % Skip to beginning of next column of original image
                index = index + 2*halfSEM;
            end
        end
        
        function out = graygradient(img, SE)

            % img: grayscale image 
            % SE:  structuring element (assumed binary)

            out = graydilate(img,SE) - grayerode(img,SE);

        end

        function Connected = findCC(img,connectivity)
            %OnePass: One-pass connected component labeling algorithm
            %   See Wikipedia page "Connected-component labeling"

            % Put a buffer of zeros around the input image
            [M,N] = size(img.Value);
            buf = zeros(M+2,N+2);
            buf(2:M+1,2:N+1) = img.Value;

            % Initialize
            [M,N] = size(buf);
            Connected = zeros(M,N);
            Mark = 1;                   % "Value"
            Difference = 1;             % "Increment"
            % Index = [];
            Nobj = 0;
            if connectivity == 4
                Offsets = [-M; -1; M; 1];
            elseif connectivity == 8
                Offsets = [-M+(-1:1) -1 1 M+(-1:1)]';
            else
                error('Connectivity should be specified as 4 or 8');
            end

            % Iterate across rows of original image pixels 
            for ii = 2:M-1
                for jj = 2:N-1
                    if buf(ii,jj)==1
                        Nobj = Nobj + 1;
                        Index = ((jj-1)*M + ii);
                        Connected(Index) = Mark;
                        while ~isempty(Index)
                            buf(Index) = 0;
                            Neighbors = bsxfun(@plus,Index,Offsets');
                            Neighbors = unique(Neighbors(:));
                            Index = Neighbors(buf(Neighbors)==1);
                            Connected(Index) = Mark;
                        end
                        Mark = Mark + Difference;
                    end
                end
            end
            Connected = Connected(2:M-1,2:N-1);
        end

    end
end


