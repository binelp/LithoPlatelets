%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ETH Zurich, Switzerland
% Separation Processes Laboratory
%
% Project:  Lithoplatelets
% Year:     2021
% MATLAB:   R2019b, Windows 64bit
% Authors:  Anna Jaeggi (AJ)
%
% Purpose:
% Simluates contours of virtual particles and tries to estimate the particles 
% orientation based on contours only.(Adapted from Shape/VTB/VTB_main.m, 
% branch: platelets, commitID: 1fca5a2)
%
% Input arguments:
% - L:                      Particle Lengths
% - realQuats:              Quaternions representing orientation of particles 
% - options:                Structure used to set options. Possible fields:
%   * crystalHeteronym:     Unique heteronym of crystal type to be modeled 
%                           (must exist in MAT file "crystalDataFileName" loaded from disk)
%   * instruments:          Cell array identifying which sizing techniques should be used/simulated
%   * nsample (opt):        Number of samples to be taken from the PSSD; {default = 2000}
%   * A (opt.):             A matrix of crystal type (fetched if not given)
%   * M (opt):              M matrix of crystal type (fetched if not given)
%   * nppimg (opt.):        Imaging methods only: number of particles per image
%   * cellDim (opt.):       Imaging methods only: dimension of the cell in pixels (Nx * Ny * Nz; z is shared axis for DP)(default: FTC)
%   * magnification (opt.): Imaging methods only: micron/pixel ratio
%   * parforProfile:        - false: use for loop for generating artificial particles
%                           - 'local': use local worker pool with options.parforNumCores
%                           - 'eulerLSF': use eulerLSF8h cluster profile with
%                           options.parforNumCores, options.parforLSFTimeLimit, and
%                           options.parforMemoryRequest
%   * parforNumCores:       Number of cores/workers to use for parfor
%   * parforLSFTimeLimit:   Time limit for the job in [h]
%   * parforMemoryRequest:  Memory requirement per core/worker in [MB]
%   * SEsizes:              Structuring element size for erosion and
%                           dilation (simulates blur)
%
% Output arguments:
% - out:                    Struct containing results, i.e. real and measured
%                           quaternions, and also options for reproducibility
% - options:                Updated options structure (-> stored A and M matrices; 
%                           reduces number of fetches from DB)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out, options] = VTB_main_rotation(L, realQuats, options)

%% Settings and input checks

% Settings and initializations
knownTechniques = {'DP'};
VTBStruct = [];
crystalDataFileName = 'CrystalData.mat';
% Check whether input is valid
options = check_InputValidity(options,knownTechniques,crystalDataFileName);
% Initialize geometry toolboxes (or download if necessary)
check_geometryToolboxes;

%% Sample PSSD and create polytopes

% uniform pssd
particleList = repmat(L,[options.nsample,1]); % [um]

% Append the list to the options
options.particleList = particleList; % [um]

% If the particles were sampled from a PSSD describing cylindrical particles, the
% dimensions (which represent distances of facets from the crystal center within the VTB!)
% need to be divided by two
if any(strcmpi(options.crystalHeteronym,{'cyl','cylinder','octarod','platelet','hexPlatelet','diagPlatelet','octaplatelet'}))
   particleList = 0.5*particleList; % [um]
end

% Create nsample polytopes using MPT3.0 and store them in an array
% If technique_IA is the only simulated measurement technique, the generation of the
% polyhedrons can actually be removed from here since they will be shifted/regenerated
% inside technique_IA.m anyway
P(options.nsample) = Polyhedron;
for i = 1:options.nsample   
    P(i) = Polyhedron(options.A,options.M*particleList(i,:).');
end
%% Simulate measurements

% Pre-allocate structure
for i = 1:length(knownTechniques) 
    VTBStruct.(knownTechniques{i}) = [];
end
% Temporarily turn of warning that occurs in 2016 only
warning('off','MATLAB:nargchk:deprecated')
% Imaging methods. Returns DP, and SP, as well as P. P is now randomly
% rotated, translated and without any particle overlap.

[out] = technique_IA_rotation(P, realQuats, options);

% Re-activate warnings
warning('on','MATLAB:nargchk:deprecated')

end

%% Auxiliary functions

function options = check_InputValidity(options,knownTechniques,crystalDataFileName)
% Check validity of input arguments to VTB_main.m

% crystalHeteronym
if ~isfield(options,'crystalHeteronym') 
    error('VTB_main:inputCheck:noCrystalHeteronym',...
        'Mandatory field ''crystalHeteronym'' was missing in options structure.');
end
% nsample
if ~isfield(options,'nsample')
    options.nsample = 2000;    % total number of particles to be simulated
    warning('VTB_main:inputCheck:noNsample',...
        ['Number of particles per measurement (ppm) was not specified. Default (',...
        num2str(options.nsample), ') option will be used.']);
end
% instruments
if ~isfield(options,'instruments')
    error('VTB_main:inputCheck:noPATsField',...
        ['Obligatory field ''instruments'' was missing in options structure.',...
        ' Please indicate what techniques should be simulated.']);
end
if ischar(options.instruments)
    options.instruments = {options.instruments};
end
if isempty(options.instruments)
    warning('VTB_main:inputCheck:nooptionsInstruments',...
        'Field options.instruments was empty, returned empty VTBStruct.');
    return
end
if any(~ismember(lower(options.instruments),lower(knownTechniques)))
    error('VTB_main:inputCheck:unknownSizingTechnique',...
        ['Unknown technique(s): ', strjoin(options.instruments(~ismember(lower(options.instruments),lower(knownTechniques))),', '),...
        '. \nThe following PAT tools are known: ',strjoin(knownTechniques,', '),'. -> Aborted.']);
end
% Assign A and M matrices if they weren't provided via the options structure
if (~isfield(options,'A') || ~isfield(options,'M')) || (isempty(options.A) || isempty(options.M))
    % Try to load crystal structure info from disk
    [crystal,~] = getCrystalInfo({options.crystalHeteronym},crystalDataFileName);
    % Assign matrices needed to construct convex polytope corresponding to the demanded
    % crystal heteronym
    options.A = crystal.H; % Matrix containing facet normal vectors in rows
    options.M = crystal.M; % Selector matrix that maps crystal dimensions to normal distances of facets to crystal center
end  
% Parallelization options
if ~isfield(options,'parforProfile')
   % Default: no parallelization -> use for-loop instead of parfor
   options.parforProfile = false;
else
   % Check if parfor profile options are valid
   validParforProfile = true;
   if islogical(options.parforProfile)
		if ~(options.parforProfile == false)
			validParforProfile = false;
		end
   else
       if ~(strcmpi(options.parforProfile,'local') || strcmpi(options.parforProfile,'eulerLSF'))
           validParforProfile = false;
       end
   end
   if ~validParforProfile
       error('VTBmain:inputCheck:invalidParforProfile',...
           'The option parforProfile must be either set to ''local'', ''eulerLSF'' or to false.');
   else
       if ~isfield(options,'parforNumCores') || ~isnumeric(options.parforNumCores)
           error('VTBmain:inputCheck:invalidParforNumCores',...
               'The option parforNumCores is required if parforProfile is set.');
       end
       if strcmpi(options.parforProfile,'eulerLSF')
           if ~isfield(options,'parforLSFTimeLimit') || ~isnumeric(options.parforLSFTimeLimit)
               error('VTBmain:inputCheck:invalidParforLSFTimeLimit',...
                   'The option parforLSFTimeLimit is required if parforProfile is ''eulerLSF''.');
           end
           if ~isfield(options,'parforMemoryRequest') || ~isnumeric(options.parforMemoryRequest)
               error('VTBmain:inputCheck:invalidParforMemoryRequest',...
                   'The option parforMemoryRequest is required if parforProfile is ''eulerLSF''.');
           end
       end
   end
end
% Plot options
if ~isfield(options,'plotResults') || islogical(options.plotResults)
    options.plotResults = false;
end

% angle Factors
if ~isfield(options,'angleFactors')
    options.angleFactors=[nan nan nan nan nan nan];
end

if ~isfield(options,'doPlots')
    options.doPlots=false;
end

if ~isfield(options,'saveFlag')
    options.saveFlag=false;
end

% Use triple projection
if ~isfield(options,'thirdCam') 
    options.thirdCam = false;
end

% Collect data for machine learning model
if ~isfield(options,'collectData') 
    options.collectData = false;
end

% plots don't appear in parallel
if options.doPlots
    options.parforProfile = false;
end

% The machine learning model only uses individual particles as datapoints
if options.collectData
    options.nppimg=1;
end

end % function


function check_geometryToolboxes
% Find or install MPT 3.0 toolbox

% Check for MPT3.0 toolbox a first time
if ~exist('mpt_init','file')
    % If MPT3.0 was not found, check for tbxmanager and run it
    if exist('tbxmanager','file')
        tbxmanager restorepath
    else
        % tbxmanager not yet found on path -> try to find it two directory levels above
        trialPathTbx = ['..',filesep','..',filesep,'tbxmanager'];
        if exist(trialPathTbx,'dir')
            addpath(trialPathTbx);
            disp('Found tbxmanager and added it to path.');
            tbxmanager restorepath
        else
            error('VTB_main:tbxmanagerNotFound',...
                ['tbxmanager.m has not been found in your search path. It is a prerequisite to run the VTB. ',...
                'Please install this file on your search path.']);
        end
    end
    % Check for MPT3.0 toolbox a second time after running the tbxmanager
    if ~exist('mpt_init','file')
        warning('VTB_main:mptNotFound',...
            'MPT3 has not been found in your search path. It and associated codes are necessary to run the eVTB.');
        usin = input('Do you want me to install all necessary toolboxes now [n/y]?','s');
        if strcmpi(usin,'y')
            tbxmanager install cddmex fourier lcp mpt
        end
    end
end

end