
function eps_to_pdf(file_path)
    status = system('export LC_CTYPE=en_US.UTF-8');
    status = system('export LC_ALL=en_US.UTF-8');

    system_command = ['find . -wholename "' file_path '" -exec epstopdf {} ";"'];
    status = system(system_command);
end
