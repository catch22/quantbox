function install_quantbox

FOLDERS = {'channels', 'info', 'linalg', 'misc', 'sdp', 'states'};
VERSION = '0.0.1-dev';

% print banner
fprintf('---------------------------------------------------------------------------\n');
fprintf('QuantBox v%s\n', VERSION);
fprintf('---------------------------------------------------------------------------\n\n');

% determine base path
basepath = fileparts(mfilename('fullpath'));
fprintf('  Directory: %s\n', basepath);

% determine environment
if exist('OCTAVE_VERSION', 'builtin')
  env = 'Octave';
else
  env = 'Matlab';
end
fprintf('  %s %s on %s\n', env, version(), computer());

% add folders to path
fprintf('\nAdding QuantBox to path:\n\n');
% addpath(basepath);
for folder = FOLDERS
  p = [basepath filesep folder{1}];
  fprintf('  %s %s \n', p, repmat('.', 1, 55 - numel(p)));
  addpath([p]);
end

fprintf('\nQuantBox has been successfully installed.\n');
fprintf('Please save the path to avoid running install_quantbox for every session.\n\n');

end
