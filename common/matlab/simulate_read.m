function [count, read] = simulate_read(haps, snp, chrom, Param)

%simulate_read    simulate different types of read for given haplotype(s)

chromInd = (snp(:, 1) == chrom);
if strcmp(Param.type, 'PE')
  readSpanLength = Param.endLen * 2 + Param.insertMean + 100 * Param.insertStd;
  readCovLength = Param.endLen * 2;
elseif strcmp(Param.type, 'Short')
  readSpanLength = Param.endLen;
  readCovLength = Param.endLen;
end

read = cell(size(haps, 2), 1);
count = cell(size(haps, 2), 1);

for hi = 1:size(haps, 2)
    
  hap = double(haps(chromInd, hi));
  hapsnp = snp(chromInd, :);    

  startIdx = hapsnp(1, 2) - readSpanLength;
  endIdx = hapsnp(end, 2) + readSpanLength;  
  
  hapIdxDict = spconvert([hapsnp(:, 2), ones(length(hap), 1), (1:length(hap))']);
  hapSnpDict = spconvert([hapsnp(:, 2), ones(length(hap), 1), hap * 2 - 1]);
  hapIdxDict(startIdx, 1) = 0;
  hapIdxDict(endIdx, 1) = 0;
  hapSnpDict(startIdx, 1) = 0;
  hapSnpDict(endIdx, 1) = 0;  
  
  nRead =  round((endIdx - startIdx) * Param.coverage / readCovLength);
  
  read{hi} = sparse(sum(chromInd), nRead);
  
  readIdx = floor(rand(1, nRead) * (endIdx - startIdx - readSpanLength)) + startIdx;
  if strcmp(Param.type, 'PE')
    insertSize = round(normrnd(Param.insertMean, Param.insertStd, 1, nRead));
    insertSize(insertSize > Param.insertMean + 100 * Param.insertStd) = ...
      Param.insertMean + 99 * Param.insertStd;
    
    for i = 1:nRead  
      readSnpIdx1 = hapIdxDict(readIdx(i) - 1 + (1:Param.endLen));
      readSnpIdx2 = hapIdxDict(readIdx(i) - 1 + Param.endLen + insertSize(i) + (1:Param.endLen));
      readSnpIdx = [readSnpIdx1; readSnpIdx2];
      if nnz(readSnpIdx) == 0; continue; end;
      
      readSnp1 = hapSnpDict(readIdx(i) - 1 + (1:Param.endLen));
      readSnp2 = hapSnpDict(readIdx(i) - 1 + Param.endLen + insertSize(i) + (1:Param.endLen));
      readSnp = [readSnp1; readSnp2];
      
      if(nnz(readSnpIdx) ~= nnz(readSnp))
        pause;
      end
      read{hi}(nonzeros(readSnpIdx), i) = nonzeros(readSnp);      
    end
  elseif strcmp(Param.type, 'Short')
    for i = 1:nRead
      readSnpIdx = hapIdxDict(readIdx(i) - 1 + (1:Param.endLen));
      if nnz(readSnpIdx) == 0; continue; end;
      
      readSnp = hapSnpDict(readIdx(i) - 1 + (1:Param.endLen));
        
      read{hi}(nonzeros(readSnpIdx), i) = nonzeros(readSnp);
    end
  end
  
  % simualte error
  nBit = nnz(read{hi});
  error = (rand(nBit, 1) < Param.error);
  read{hi}(error) = -read{hi}(error);
  
  % read count for each snp
  count{hi} = [sum(read{hi} == -1, 2), sum(read{hi} == 1, 2)];
end