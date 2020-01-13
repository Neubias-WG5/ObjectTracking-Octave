#! /bin/octave -qf

SharePenalty = 0.5; 
AlphaCM = 0.33;
SmallPart = 15;
OvlItExtend = 2;
MaxRescueDistance = 48;

arg_list = argv ();
InputFolder = arg_list{1};
OutputFolder = arg_list{2};
GBlurRad = str2num(arg_list{3}); 
MeanBox = str2num(arg_list{4});
Offset = str2num(arg_list{5});
if nargin > 5
	SharePenalty = str2num(arg_list{6}); 		% Penalty for sharing event (split particle) in likelihood metric (higher --> less likely to happen), default: 0.5
	AlphaCM = str2num(arg_list{7});		 		% Weight for centre of mass displacement in likelihood metric (higher --> more weigth), default: 0.33
	SmallPart = str2num(arg_list{8});	 		% Remove particles of area smaller than SmallPart, set to Inf to disable, default: 15
	OvlItExtend = str2num(arg_list{9});	 		% Expand object edges by this amount to allow further overlap (experimental), default: 2
	MaxRescueDistance = str2num(arg_list{10});	% Maximum distance to rescue object, set to -1 to prevent object rescue
end

%% Load image toolbox
pkg load image; 

se = [0 1 0;1 1 1;0 1 0];    % Structure element used for overlap map dilations (outside objects and inside particles)
se2 = [0 1 0;1 1 1;0 1 0];   % Structure element used for object map closing (merge almost touching objects with identical labels) 

UnassignState = 0;          % Add a new particle state: particle not assigned to any object (experimental)
Reassign = 0;               % Try to assign new object index to a free index
Debug = 0;                  % Display debg information

Graph = 0;          % Plot results 
Export = 1;         % Export result images to files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Files = dir(strcat([InputFolder '*.tif']));
num_movies = numel(Files);

for m = 1:num_movies

inf = imfinfo(strcat([InputFolder Files(m).name]));
num_images = numel(inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Graph == 1
    FigParticles = figure;
    FigOverlap = figure;
    FigObjects = figure;
    Yplot = ceil(sqrt(num_images));
    Xplot = ceil((num_images)/Yplot);
end
AccBestDiffArea = 0;
AccBestDiffCM = 0;
DivLog = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ObjLblMask = uint16(zeros(inf(1).Height,inf(1).Width,1,num_images));

for kf = 1:num_images

    fprintf('.');
    
	%% Read image
	A = imread(strcat([InputFolder Files(m).name]),kf);
	
    %% Threshold image (+ split objects)
	A = imsmooth(single(A),'Gaussian',GBlurRad,GBlurRad) >= (Offset+imfilter(single(A),ones(MeanBox,MeanBox)/(MeanBox*MeanBox),'same','symmetric'));
	if kf > 1
	  D = -bwdist(~A);
	  D = imsmooth(D,'Gaussian',1,1);
	  Ld = watershed(D);
	  A(Ld == 0) = 0;
	end
	if ~(SmallPart == Inf)
		A = bwareaopen(A, SmallPart)*255;
	end
	%A = imclearborder(A);   % Remove objects touching edges

    %% Connected particles at current iteration 
    CC = bwconncomp(A == 255, 8);
    Particles = CC.PixelIdxList;    
    Npart = length(Particles);

    %% Fill particle label map + compute particle statistics
    PartLbl = zeros(size(A));
    AreaPart = zeros(1,Npart);
    CMPartX = zeros(1,Npart); 
    CMPartY = zeros(1,Npart);
    for i = 1:Npart
        PartLbl(Particles{i}) = i; 
        AreaPart(i) = length(Particles{i});
        [Y X] = ind2sub(size(A),Particles{i});
        CMPartX(i) = mean(X);
        CMPartY(i) = mean(Y); 
    end

    if kf>1

        %% Compute zone of influence of the objects outside particles (masked dilation)   
        if OvlItExtend>0
            ObjLblExt = ObjLbl;   
            for i = 1:OvlItExtend
                ObjLblExt2 = imdilate(ObjLblExt,se);
                ObjLblExt = (ObjLblExt==0).*ObjLblExt2+(ObjLblExt>0).*ObjLblExt;
            end
        end

        %% Compute particle overlap map
        OvlMap = zeros(size(A));
        if OvlItExtend>0
            for i = 1:Npart
                OvlMap(Particles{i}) = ObjLblExt(Particles{i});
            end
        else
            for i = 1:Npart
                OvlMap(Particles{i}) = ObjLbl(Particles{i});
            end
        end

        %% Expand zone of influence of the objects inside particles (until convergence)
        test = 1;
        %cnt = 0;
        OvlMapOld = OvlMap;
        while test>0
            OvlMapBuf = imdilate(OvlMap,se).*(PartLbl > 0);
            OvlMapNew = (OvlMap == 0).*OvlMapBuf+(OvlMap > 0).*OvlMap;
            Dif = OvlMapNew-OvlMap;
            test = any(Dif(:));
            %cnt = cnt+1;
            OvlMap = OvlMapNew;
        end 

    end

    if (kf == 1)
        %% All particles --> objects (first frame)
        ObjLbl = PartLbl;
        AreaObj = AreaPart;
        LastAreaObj = AreaObj;
        CMObjX = CMPartX;
        CMObjY = CMPartY;
        LastCMObjX = CMObjX;
        LastCMObjY = CMObjY;
        Nobj = Npart;
        ScaleImg = Nobj;

        %% Display object label map
        if Graph==1
            figure(FigObjects);
            subplot(Xplot,Yplot,1);
            imagesc(ObjLbl,[0 ScaleImg]);
            set(gcf,'name',strcat('Initial Objects ',num2str(kf)));
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
        end

    else

        %% Display label maps
        if ((Graph == 1) && (FigParticles > -1))
            figure(FigParticles);
            imagesc(PartLbl,[0 ScaleImg]);
            set(gcf,'name',strcat('Particles ',num2str(kf)));
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
        end

        %% Export overlap and particle label maps
        if ((Graph == 1) && (FigOverlap > -1))
            figure(FigOverlap);
            imagesc(OvlMap,[0 ScaleImg]);
            set(gcf,'name',strcat('OvlMap ',num2str(kf)));
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
        end    

        %% Find subparticles (same object) and compute their statistics 
        PartObjIndx = cell(1,Npart);
        PartObjArea = cell(1,Npart);
        PartObjCMX = cell(1,Npart);
        PartObjCMY = cell(1,Npart);
        for i = 1:Npart
            Indx = OvlMap(Particles{i}).';
            UniqueIndx = unique(Indx);
            PartObjIndx{i} = UniqueIndx;
            AreaVector = zeros(1,length(UniqueIndx));
            CMXVector = zeros(1,length(UniqueIndx));
            CMYVector = zeros(1,length(UniqueIndx));
            for k =1:length(UniqueIndx)
                Msk = (Indx == UniqueIndx(k));
                AreaVector(k) = sum(Msk);
                [Y X] = ind2sub(size(A),Particles{i}(Msk));
                CMXVector(k) = mean(X);
                CMYVector(k) = mean(Y);
            end
            PartObjArea{i} = AreaVector;
            PartObjCMX{i} = CMXVector;
            PartObjCMY{i} = CMYVector;
        end

        %% For each object: Find index of overlapping particle(s)
        ObjPartIndx = cell(1,Nobj);
        for i = 1:Npart
            Obj = PartObjIndx{i};
            for j = 1:length(Obj)
                if (Obj(j) > 0)
                    ObjPartIndx{Obj(j)} = [ObjPartIndx{Obj(j)} i];
                end
            end
        end
        for i = 1:Nobj
            ObjPartIndx{i} = unique(ObjPartIndx{i});
        end

        %% Group object(s) in clusters so that particle(s) of new frame only overlap objects of a given cluster
        %% Important to speed up the optimization
        ExploredObj = zeros(1,Nobj); % Binary flag to state that the object has already been considered
        ObjCluster = cell(1,Nobj);   % Object(s) making up each cluster 
        PartCluster = cell(1,Nobj);  % Particle(s) overlapping at least one object of a given cluster
        NCluster = 0;
        for i = 1:Nobj
            if ((length(ObjPartIndx{i}) > 0)&&(ExploredObj(i) == 0)) 
                NCluster = NCluster+1;
                CurrentObj = i;
                CurrentPart = [];
                NewObj = i;
                test = 0;
                while (test == 0)
                    WanderedObj = [];
                    for j = 1:length(NewObj)
                        PartInvolved = ObjPartIndx{NewObj(j)};
                        for k = 1:length(PartInvolved)
                            WanderedObj = [WanderedObj PartObjIndx{PartInvolved(k)}];
                        end
                        WanderedObj = unique(WanderedObj);
                        CurrentPart = unique([CurrentPart PartInvolved]);
                    end
                    NewObj = setdiff(WanderedObj,CurrentObj);
                    test = isempty(NewObj);
                    CurrentObj = unique([CurrentObj WanderedObj]);
                end
                ObjCluster{NCluster} = CurrentObj;
                PartCluster{NCluster} = CurrentPart;
                ExploredObj(CurrentObj) = 1;
            end
        end

        %% Optimize object/particle association within each cluster 
        %% NOTE: Possible improvement --> include a state "no object wins the particle"
        NewAreaObj = zeros(1,Nobj);
        NewCMXObj = zeros(1,Nobj);
        NewCMYObj = zeros(1,Nobj);
        NewMultObj = zeros(1,Nobj);
        BestObjArea = zeros(1,Nobj);
        BestObjCMX = zeros(1,Nobj);
        BestObjCMY = zeros(1,Nobj);
        BestMultObj = zeros(1,Nobj);
        PartAssign = zeros(1,Npart);
        BestPartAssign = zeros(1,Npart);

        for i = 1:NCluster

            %% Cluster optimization initialization
            BestDiff = Inf;
            BestState = -1;
            ObjSet = ObjCluster{i};
            PartSet = PartCluster{i};
            NPartSet = length(PartSet);
            TotPartStates = zeros(1,NPartSet);
            for j = 1:length(PartSet)
                L = length(PartObjIndx{PartSet(j)}); % L is the number of object(s) the particle overlaps with
                TotPartStates(j) = (L+(L>1)+1*(UnassignState == 1)); 
                % First L states: One object wins the particle (L first states since L objects in cluster) 
                % (L+1)th state: All objects share the particle 
                % (L+2)th state: No object wins the particle (optional)
            end

            %% Cluster optimization: test all possible assignements
            % PartStates: A matrix with as many line as possible associations, one column per particle, encoding particle state
            NCnt = length(TotPartStates);
            NComb = prod(TotPartStates);
            PartStates = zeros(NComb,NCnt);
            Cnt = ones(1,NCnt);
            Cnt(NCnt) = Cnt(NCnt)-1;
            for k = 1:NComb 
                Cnt(NCnt) = Cnt(NCnt)+1;
                for z = NCnt:-1:1
                    if (Cnt(z) > TotPartStates(z))
                    Cnt(z) = 1;
                    Cnt(z-1) = Cnt(z-1)+1;
                    end
                end
            PartStates(k,:) = Cnt;
            end

            %% Compute likelihood metric for all possible particle to object associations
            for j = 1:size(PartStates,1)   
                CurrentState = PartStates(j,:);     
                NewAreaObj(ObjSet) = 0;
                NewCMXObj(ObjSet) = 0;
                NewCMYObj(ObjSet) = 0;
                NewMultObj(ObjSet) = 0;
                PartAssign(PartSet) = 0;
                StatePenalty = 0;
                for k = 1:NPartSet
                    ObjSubSetIndx = PartObjIndx{PartSet(k)};
                    ObjSubSetArea = PartObjArea{PartSet(k)};
                    ObjSubSetCMX = PartObjCMX{PartSet(k)};
                    ObjSubSetCMY = PartObjCMY{PartSet(k)};
                    NObjSubSet = length(ObjSubSetIndx);
                    if (CurrentState(k) == TotPartStates(k)) % Share particle between different objects
                        PartAssign(PartSet(k)) = 0;
                        if (NObjSubSet == 1) % Single object inside particle
                            NewAreaObj(ObjSubSetIndx(1)) = NewAreaObj(ObjSubSetIndx(1))+AreaPart(PartSet(k));
                            NewCMXObj(ObjSubSetIndx(1)) = CMPartX(PartSet(k));
                            NewCMYObj(ObjSubSetIndx(1)) = CMPartY(PartSet(k));
                            NewMultObj(ObjSubSetIndx(1)) = NewMultObj(ObjSubSetIndx(1))+1;
                            PartAssign(PartSet(k)) = ObjSubSetIndx(1);
                        else
                            for l = 1:NObjSubSet % General case: each object takes its own share
                                NewCMXObj(ObjSubSetIndx(l)) = (NewCMXObj(ObjSubSetIndx(l))*NewAreaObj(ObjSubSetIndx(l))+ObjSubSetCMX(l)*ObjSubSetArea(l))/(NewAreaObj(ObjSubSetIndx(l))+ObjSubSetArea(l));
                                NewCMYObj(ObjSubSetIndx(l)) = (NewCMYObj(ObjSubSetIndx(l))*NewAreaObj(ObjSubSetIndx(l))+ObjSubSetCMY(l)*ObjSubSetArea(l))/(NewAreaObj(ObjSubSetIndx(l))+ObjSubSetArea(l));
                                NewAreaObj(ObjSubSetIndx(l)) = NewAreaObj(ObjSubSetIndx(l))+ObjSubSetArea(l);
                                NewMultObj(ObjSubSetIndx(l)) = NewMultObj(ObjSubSetIndx(l))+1;
                            end
                            StatePenalty = StatePenalty + SharePenalty*AreaPart(PartSet(k));
                        end         
                    else  
                        if ( ( CurrentState(k) == (TotPartStates(k)-1) )&&( UnassignState == 1 ) )
                            % No object wins the particle
                            PartAssign(PartSet(k)) = NaN;
                        else
                            % One object wins the particle
                            PartAssign(PartSet(k)) = ObjSubSetIndx(CurrentState(k));
                            NewAreaObj(ObjSubSetIndx(CurrentState(k))) = NewAreaObj(ObjSubSetIndx(CurrentState(k)))+AreaPart(PartSet(k));
                            NewCMXObj(ObjSubSetIndx(CurrentState(k))) = CMPartX(PartSet(k));
                            NewCMYObj(ObjSubSetIndx(CurrentState(k))) = CMPartY(PartSet(k));
                            NewMultObj(ObjSubSetIndx(CurrentState(k))) = NewMultObj(ObjSubSetIndx(CurrentState(k)))+1;
                        end
                    end
                end

                %% Compute likelihood metric for this state (smallest metric --> best)
                DiffArea = sum(abs(AreaObj(ObjSet)-NewAreaObj(ObjSet)));
                DiffCM = sum(((CMObjX(ObjSet)-NewCMXObj(ObjSet)).*(NewAreaObj(ObjSet)>0)).^2+((CMObjY(ObjSet)-NewCMYObj(ObjSet)).*(NewAreaObj(ObjSet)>0)).^2);
                Diff = DiffArea+AlphaCM*DiffCM+StatePenalty;
                if Diff<BestDiff
                    BestState = CurrentState;
                    BestDiff = Diff;
                    BestDiffArea = DiffArea;
                    BestDiffCM = DiffCM;
                    BestPartAssign(PartSet) = PartAssign(PartSet);
                    BestObjArea(ObjSet) = NewAreaObj(ObjSet);
                    BestObjCMX(ObjSet) = NewCMXObj(ObjSet);
                    BestObjCMY(ObjSet) = NewCMYObj(ObjSet);
                    BestMultObj(ObjSet) = NewMultObj(ObjSet);
                end
            end

            if (Debug == 1)
                disp(strcat('Cluster_',num2str(i),'_NStates_',num2str(prod(TotPartStates)),'____PartIndx_',num2str(PartSet),'_PartArea_',num2str(AreaPart(PartSet)),'____ObjIndx_',num2str(ObjSet),'_ObjArea_',num2str(AreaObj(ObjSet)))) ;
                StrState = 'Best State: ';
                for z = 1:length(BestState)
                    StrState = strcat(StrState,num2str(BestState(z)));
                end
                disp(StrState);
                disp(strcat('AssignedObjects_',num2str(BestPartAssign(PartSet)),'_AreaObjects_',num2str(BestObjArea(ObjSet)),'_Metric_',num2str(BestDiffArea),'_',num2str(round(BestDiffCM)),'_',num2str(BestDiff)));
            end

            % Combined metric
            AccBestDiffCM = AccBestDiffCM+BestDiffCM;
            AccBestDiffArea = AccBestDiffArea+BestDiffArea;
        end

        %% Fill object label map
        AreaObj = BestObjArea;
        CMObjX = BestObjCMX;
        CMObjY = BestObjCMY;
        ObjLbl = OvlMap;
        for i = 1:Npart
            if(BestPartAssign(i)>0)
                ObjLbl(Particles{i}) = BestPartAssign(i);
            end
            if(BestPartAssign(i) == NaN)
                ObjLbl(Particles{i}) = 0;
            end
        end

        %% Detect object disappearances
        Disappeared = (((LastAreaObj > 0)-(AreaObj > 0)) > 0);
        DisappearedIndx = find(Disappeared == 1);

        %% Handle object rescues and appearances
        DiffMask = (PartLbl>0)-(ObjLbl>0);
        if any(any(DiffMask,1))
            CC = bwconncomp(DiffMask, 8);
            Rescued = zeros(1,CC.NumObjects);
            NewArea = zeros(1,CC.NumObjects);
            NewCMObjX = zeros(1,CC.NumObjects);
            NewCMObjY = zeros(1,CC.NumObjects);
            for i=1:CC.NumObjects
                NewArea(i) = length(CC.PixelIdxList{i});
                [Y X] = ind2sub(size(A),CC.PixelIdxList{i});
                NewCMObjX(i) = mean(X);
                NewCMObjY(i) = mean(Y);
            end

            %% Handle object rescues  
            %% Note: The assignement could be optimized
            for i=1:CC.NumObjects 
                IndxObj = 0;    
                MinDist = Inf;
                Minj = [];
                if any(Disappeared)
                    for j=1:length(DisappearedIndx)
                        Dist = sqrt((LastCMObjX(DisappearedIndx(j))-NewCMObjX(i))^2+(LastCMObjY(DisappearedIndx(j))-NewCMObjY(i))^2);
                        if ((Dist < MaxRescueDistance)&&(Dist < MinDist))
                            IndxObj = DisappearedIndx(j);
                            MinDist = Dist;
                            Minj = j;
                        end
                     end
                end
                if (IndxObj > 0) % Rescued
                    ObjLbl(CC.PixelIdxList{i}) = IndxObj;
                    AreaObj(IndxObj) = NewArea(i);
                    CMObjX(IndxObj) = NewCMObjX(i);
                    CMObjY(IndxObj) = NewCMObjY(i);
                    BestMultObj(IndxObj) = 1;
                    Rescued(i) = 1;
                    DisappearedIndx(Minj)=[];
                end
            end

            %% Handle object appearances
            for i=1:CC.NumObjects
                if (Rescued(i) == 0)
                    if (Reassign == 1) % Try to assign new object index to a free index
                        Candidates = find(AreaObj==0); % Find empty object index
                        if isempty(Candidates) % No empty object index --> increase object counter
                            Nobj = Nobj+1;
                            IndxObj = Nobj;
                        else
                            IndxObj = Candidates(1);
                        end
                    else % Always increase object counter
                        Nobj = Nobj+1;
                        IndxObj = Nobj;
                    end
                    ObjLbl(CC.PixelIdxList{i}) = IndxObj;
                    AreaObj(IndxObj) = NewArea(i);
                    CMObjX(IndxObj) = NewCMObjX(i);
                    CMObjY(IndxObj) = NewCMObjY(i);
                    BestMultObj(IndxObj) = 1;
                end
            end
        end

        %% Handle object divisions and mergings
        if ~(isempty(find(BestMultObj>1)))
            ObjLblClosed = imclose(ObjLbl,se2);
            for i=1:Nobj
                if BestMultObj(i)>1 % Various particles assigned to the same object  
                    CC = bwconncomp((ObjLblClosed == i), 8); 
                    if (CC.NumObjects > 1) % Particles are not closely touching: candidate division
                        %% Find the closest particle to the original object
                        MinDist = Inf;
                        MinDistIndx = 0;
                        for j=1:CC.NumObjects
                            [Y X] = ind2sub(size(A),CC.PixelIdxList{j});
                            Dist = sqrt((mean(X)-LastCMObjX(i))^2+(mean(Y)-LastCMObjY(i))^2);
                            if(Dist<MinDist)
                                MinDistIndx = j;
                                MinDist = Dist;
                            end
                        end
                        for j=1:CC.NumObjects
                            %% Update area and position of mother object
                            if (j == MinDistIndx)
                                %IndxObj = i;                       %% Modified so that both children have new indices
                                Nobj = Nobj+1;                      %%
                                IndxObj = Nobj;                     %%
                                AreaObj(IndxObj) = length(CC.PixelIdxList{j});
                                [Y X] = ind2sub(size(A),CC.PixelIdxList{j});
                                CMObjX(IndxObj) = mean(X);
                                CMObjY(IndxObj) = mean(Y);
                                %ObjLbl(CC.PixelIdxList{j}) = i;    %% Modified so that both children have new indices
                                ObjLbl(CC.PixelIdxList{j}) = IndxObj;
                                Div1 = IndxObj;
                            else
                            %% Division --> new object (not same index to avoid re-merge)    
                                AreaCand = length(CC.PixelIdxList{j});
                                Nobj = Nobj+1;
                                IndxObj = Nobj;
                                ObjLbl(CC.PixelIdxList{j}) = IndxObj;
                                AreaObj(IndxObj) = AreaCand;
                                [Y X] = ind2sub(size(A),CC.PixelIdxList{j});
                                CMObjX(IndxObj) = mean(X);
                                CMObjY(IndxObj) = mean(Y);
                                Div2 = IndxObj;
                            end
                        end
                        disp(sprintf('In frame %s Obj %s gives birth to Obj %s and Obj %s',num2str(kf),num2str(i),num2str(Div1),num2str(Div2)));
                        %DivLog = [DivLog ; [kf i Div1 Div2]];
						Parent(Div1) = i;
						Parent(Div2) = i;
                    %% Candidates are closely touching: candidate merging        
                    else 
                         ObjLbl(CC.PixelIdxList{1}) = i;
                         AreaObj(i) = length(CC.PixelIdxList{1});
                         [Y X] = ind2sub(size(A),CC.PixelIdxList{1});
                         CMObjX(i) = mean(X);
                         CMObjY(i) = mean(Y);
                    end
                end
            end
        end

        %% Buffer current states
        LastAreaObj = AreaObj;
        LastCMObjX = CMObjX;
        LastCMObjY = CMObjY;

        %% Display object label map
        if Graph==1
            figure(FigObjects);
            subplot(Xplot,Yplot,kf-startframe+1);
            imagesc(ObjLbl,[0 ScaleImg]);
            set(gcf,'name',strcat('Objects ',num2str(kf)));
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
        end

    end

    %% Store results to object label mask
    ObjLblMask(:,:,1,kf) = uint16(ObjLbl);
    
end

if Export == 1
	imwrite(ObjLblMask,strcat(OutputFolder,Files(m).name));
end

%% Create CTC division text file
tracklives = uint8(zeros(Nobj,65536));
for i = 1:num_images
	tracklives(i,:) = imhist(ObjLblMask(:,:,1,i),65536)>0;
end
trackstart = zeros(Nobj,1);
trackend = zeros(Nobj,1);
Parent(Nobj+1) = 0;
for i = 1:Nobj
	inds = find(tracklives(:,i+1)>0);
	trackstart(i) = inds(1);
	trackend(i) = inds(end);
end
M = [(1:Nobj).' trackstart trackend Parent(1:end-1).'];

%% Write division text file
FileName = Files(m).name;
dlmwrite(strcat(OutputFolder,FileName(1:end-4),'.txt'),M,'delimiter',' ');
disp(' ');

end