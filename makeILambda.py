import os
import shutil

def make_ILambda(input_folder, file_name, destination_folder, num_copies):
    # Construct the full path of the input file
    input_file = os.path.join(input_folder, file_name)
    
    # Check if the input file exists
    if not os.path.isfile(input_file):
        print(f"File '{input_file}' not found!")
        return
    
    # Ensure the destination folder exists
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
    
    for i in range(0, num_copies ):
        # Generate a new file name with the format ILambda_{i}_0 (no extension)
        new_file_name = f"ILambda_{i}_0"
        destination_file = os.path.join(destination_folder, new_file_name)
        
        # Check if the file already exists and remove it
        if os.path.exists(destination_file):
            os.remove(destination_file)
            print(f"Deleted existing file: {destination_file}")
        
        # Copy the file
        shutil.copy(input_file, destination_file)
        print(f"Created copy: {destination_file}")


# Example usage
input_file_zoneone = 'case5_Steady_rad_sir2/0/zoneone'
destination_folder_zoneone = 'case5_Steady_rad_sir2/0/zoneone'
input_file_zonetwo = 'case5_Steady_rad_sir2/0/zonetwo'
destination_folder_zonetwo = 'case5_Steady_rad_sir2/0/zonetwo'
num_copies_zoneone = 64
num_copies_zonetwo = 64
file_name = 'ILambda'

make_ILambda(input_file_zoneone, file_name, destination_folder_zoneone, num_copies_zoneone)

make_ILambda(input_file_zonetwo, file_name, destination_folder_zonetwo, num_copies_zonetwo)
