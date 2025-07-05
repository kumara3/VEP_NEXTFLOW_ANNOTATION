#icav2 projects create test_v1 --billing-mode "Tenant" --short-descr "Test creation of pipeline for running pipelines"
#icav2 projectdata delete /test_v1_samples/
#icav2 projectdata delete /test_v1/

#!/bin/bash
check_ica_cli() {
    if ! command -v icav2 &> /dev/null; then
        echo "Error: ICA CLI not found. Please install and configure the ICA CLI first." >&2
	echo "Installation instructions:https://help.ica.illumina.com/command-line/latest/install" >&2
	return 1
    fi
    return 0
}

create_project(){
	local project_name=$1
	echo "Searching ICA..."
	project_id=$(icav2 projects list|grep -w "$project_name"|cut -f2)
	echo "$project_id"
	if [ -n "$project_id" ]; then
		echo "$project_id"
		echo "Project $project_name already exists. Skipping creation."
		return 0
	#else
        # 	echo "Project "$project_name" is not found."
	#	return 1
        fi

	echo "Creating new ICA project:$project_name..."
        if icav2 projects create $project_name --billing-mode "Tenant" --short-descr "WGS pipeline Test"; then 
	    	echo "Project $project_name create successfully"
	else
		echo "Error: Failed to create project "$project_name"" >&2
		return 1
	fi
}

upload_folder() {
    local source_path=$1
    local dest_name=$2
    dest_name_new=$(echo "$dest_name"|sed 's/\///g')
    echo $dest_name_new
    
    
    if [ ! -e "$source_path" ]; then
        echo "Error: Path not found: $source_path" >&2
        return 1
    fi
    
    echo "Uploading folder '$source_path' to ICA ..."

    project_id=$(icav2 projects list|grep -w "$dest_name_new"|cut -f1)
    echo "$project_id"
   
    if [ -z "$project_id" ]; then
	echo "Error: Project id "$project_id" for "$dest_name_new" not found in ICA." >&2
	return 1
    fi

    if icav2 projectdata upload "$source_path" "$dest_name" --project-id "$project_id";then
	echo "Successfully uploaded folder to ICA "$dest_name""
	echo "Upload loaction:$det_name in ICA"
	return 0
     else
	echo "Error: Upload failed" >&2
	return 1
    fi
    
}

main() {
    if [ $# -lt 3 ]; then
        echo "Usage: $0 <project_name><source_path><dest_path_in_ICA>" >&2
        exit 1
    fi
    
    local project_name=$1
    local source_path=$2
    local dest_name=$3

    
    if ! check_ica_cli; then
        exit 1
    fi
    
    if ! create_project "$project_name"; then
	exit 1
    fi
    
    upload_folder "$source_path" "$dest_name"
    exit $?
}

main "$@"
