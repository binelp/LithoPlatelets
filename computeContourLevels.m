function [x,y] = computeContourLevels(PSSD)

[L1, L2, L3] = PSSD.PSSD.grid.y;
q_cell = {PSSD.PSSD.F};

fullmom = compute3DCrossMoments(PSSD.PSSD);

if ~iscell(L1) || ~iscell(L2) || ~iscell(L3)
    L1Vec = L1;
    L2Vec = L2;
    L3Vec = L3;
    L1 = cell(1,length(q_cell));
    L2 = cell(1,length(q_cell));
    L3 = cell(1,length(q_cell));
    for ii=1:length(q_cell)
        L1{ii} = [L1Vec(1)-1E5*eps L1Vec L1Vec(end)+1E5*eps];
        L2{ii} = [L2Vec(1)-1E5*eps L2Vec L2Vec(end)+1E5*eps];
        L3{ii} = [L3Vec(1)-1E5*eps L3Vec L3Vec(end)+1E5*eps];
    end
end

temp_Cell   = cell(1,1);
temp_Cell{1}= permute(q_cell{1},[2 1 3]);
q_cell      = temp_Cell;


target = [0.8 0.5 0.25];
guess = [0.1 0.2 0.5];
for kk=1:length(target)
    x(kk) = fzero(@(lvl)calcRatio(lvl,{q_cell,ii,PSSD,fullmom,target(kk)}),guess(kk));
end


 for ignoredDim=1:3
        dims = find([1,2,3]~=ignoredDim);
        q_m  = {squeeze(sum(PSSD.PSSD.F,ignoredDim))};
        
        PSSD_2D.PSSD.F = q_m{1};
        PSSD_2D.PSSD.grid = PSSD.PSSD.grid(dims); %#ok<FNDSB>
        fullmom = computeCrossMoments(PSSD_2D.PSSD);
        
        for kk=1:length(target)
            y{ignoredDim}(kk) = fzero(@(lvl)calcRatio2(lvl,{q_m,ii,PSSD_2D,fullmom,target(kk)}),x(kk));
        end
        
 end

end

function R = calcRatio(lvl,p)

q_cell = p{1};
ii = p{2};
PSSD = p{3};
fullmom = p{4};
target = p{5};

normqcell = q_cell{ii}./max(q_cell{ii}(:));
logicalq = normqcell>lvl;
subPSSD.F = q_cell{1};
subPSSD.F(~logicalq) = 0;
subPSSD.grid = PSSD.PSSD.grid;
subPSSD.F = permute(subPSSD.F,[2 1 3]);
submom = compute3DCrossMoments(subPSSD);
R = submom.mu000/fullmom.mu000-target;
end

function R = calcRatio2(lvl,p)

q_cell = p{1};
ii = p{2};
PSSD = p{3};
fullmom = p{4};
target = p{5};

normqcell = q_cell{ii}./max(q_cell{ii}(:));
logicalq = normqcell>lvl;
subPSSD.F = q_cell{1};
subPSSD.F(~logicalq) = 0;
subPSSD.grid = PSSD.PSSD.grid;
submom = computeCrossMoments(subPSSD);
R = submom.mu00/fullmom.mu00-target;
end
