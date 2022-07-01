function [quatList] = sampleLithoplateletQD(nsample, L, QDpath)
    % Interpolate between experimentally measured quaternion distributions
    % (QD)
    load(QDpath)
    [QDmirrored] = interpolateQD(allQDs, L);
    
    % Undo mirroring 
    QD.F = zeros([20, 20, 20]);
    QD.F(:,11:end,:) = QDmirrored.F;
    QD.F(:,1:10,:) = flip(flip(QDmirrored.F(:,:,:),2),1);

    QD.grid = QDmirrored.grid;
    QD.grid(2).y = QD.grid(1).y ;
    QD.grid(2).boundaries = QD.grid(1).boundaries;
    
    % Sample obtained distributions
    [quatList3D] = PSSD2ParticleList(QD,nsample);
    
    % Convert reduced 3D quaternions back to 4D
    quatList = quat3Dto4D(quatList3D);
end
function [QDinterp] = interpolateQD(allQDs, L)
AS = cat(1,allQDs.AS);
QDs = {allQDs.QD};
% New point with unknown angle distribution (aspect ratio)
P = [L(1)/L(2),L(2)/L(3)];

ixs={};
% Separate known points into quadrants
ixs{2,2} = find(AS(:,1) >= P(1) & AS(:,2) >= P(2) ); % upper right
ixs{1,1} = find(AS(:,1) <= P(1) & AS(:,2) <= P(2) ); % lower left
ixs{1,2} = find(AS(:,1) <= P(1) & AS(:,2) >= P(2) ); % upper left
ixs{2,1} = find(AS(:,1) >= P(1) & AS(:,2) <= P(2) ); % lower right
% Only closest in each quadrant
for ii = 1:2
  for jj = 1:2  
    [~,minIx]=min(vecnorm(P-AS(ixs{ii,jj},:),2,2));
    ixs{ii,jj} = ixs{ii,jj}(minIx);
  end
end
% If exactly in between two pops
if length(unique([ixs{:}]))==2
    if ixs{1,2} ~= ixs{1,1} 
        denom = AS(ixs{1,2},2) - AS(ixs{1,1},2);
        I = (AS(ixs{1,2},2)-P(2))/denom*QDs{ixs{1,1}}.F  + (P(2)-AS(ixs{1,1},2))/denom*QDs{ixs{1,2}}.F ;
    else
        denom = AS(ixs{2,1},1) - AS(ixs{1,1},1);
        I = (AS(ixs{2,1},1)-P(1))/denom*QDs{ixs{1,1}}.F  + (P(1)-AS(ixs{1,1},1))/denom*QDs{ixs{2,1}}.F ;
    end
% If exactly on one pop
elseif length(unique([ixs{:}]))==1
        I = QDs{ixs{1,1}}.F ;
% General case, between 4 pops
else
    denom = AS(ixs{2,2},1) - AS(ixs{1,2},1);
    I_1 = (AS(ixs{2,2},1) - P(1))/(denom)*QDs{ixs{1,2}}.F + (P(1)-AS(ixs{1,2},1))/(denom)*QDs{ixs{2,2}}.F;
    I_0 = (AS(ixs{2,2},1) - P(1))/(denom)*QDs{ixs{1,1}}.F + (P(1)-AS(ixs{1,2},1))/(denom)*QDs{ixs{2,1}}.F;

    denom = AS(ixs{1,2},2) - AS(ixs{1,1},2);
    I = (AS(ixs{1,2},2)-P(2))/denom*I_0 + (P(2)-AS(ixs{1,1},2))/denom*I_1;
end

QDinterp.F = I;
QDinterp.grid = QDs{1}.grid;
end