function eqInitPop = InitPopulation(N,bnd)

  D       = size(bnd,1); % Read the number of dimensions
  EpD     = round(N^(1/D)); % Calulate the elements per dimension

  % Some variables
  dims    = 1 : D;  s0vars  = '';
  for iD  = 1 : D,
      s0vars  = [s0vars, sprintf('linspace(bnd(%d,1),bnd(%d,2),EpD),',iD,iD)];
  endfor
  s1vars  = sprintf('X%d,',dims); % String of elements
  s2vars  = sprintf('X%d(:),',dims); % String of elements

  eval(['[',s1vars(1:end-1),']=ndgrid(',s0vars(1:end-1),');']);

  % Get the initial population
  InitPop = eval(['[',s2vars(1:end-1),']']);
