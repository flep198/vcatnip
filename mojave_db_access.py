import mysql.connector as sql
import pandas as pd
import numpy as np
import requests
import os
import pexpect
from graph_generator import getComponentInfo

def upload_csv_to_MOJAVE(csv_file,observer,password,source):
    df=pd.read_csv(csv_file)

    #connect to database:
    mydb = sql.connect(
        host="mojavedb.mpifr-bonn.mpg.de",
        user="agn",
        password=password,
        database="galaxies"
    )

    cursor = mydb.cursor()

    for index, row in df.iterrows():
        freq=float(row["freq"])/1e9
        if 15<freq and freq<16:
            database="components"
        elif 22< freq and freq<25:
            database="Kcomponents"

        elif 42< freq and freq<45:
            database="Qcomponents"
        else:
            database="components"
        freq = "{:.1f}".format(freq)

        #check if a modelfit already exists:
        select_query = ("SELECT COUNT(*) FROM " + database + " WHERE source='" + source + "' AND observer='" + observer +
                        "' AND method='UV' AND stokes='I' AND epoch='" + row['date'] +"' AND id='" + str(row["component_number"]) + "'" )

        cursor.execute(select_query)
        result = cursor.fetchone()

        if result[0] == 0:
            # No matching row found, so insert new record
            insert_query = ("INSERT INTO " + database +
            " (source, epoch, id, flux, dist, pa, size, ratio, cpa, stokes,observer,method,freq, use_in_fit, counterjet, rating) VALUES " +
            "(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)")
            cursor.execute(insert_query, (source, row["date"], row['component_number'], row["flux"],row["radius"], row["theta"],
                                          row["size"],row["ratio"],"0","I",observer,"UV",freq,"1","-1","1"))
        else:
            update_query = ("UPDATE " + database + " SET source = %s, " +
                                                   "epoch = %s, " +
                                                   "id = %s, " +
                                                   "flux= %s, " +
                                                   "dist= %s, " +
                                                   "pa= %s, " +
                                                   "size= %s, " +
                                                   "ratio= %s, " +
                                                   "cpa= %s, " +
                                                   "stokes= %s, " +
                                                   "observer= %s, " +
                                                   "method= %s, " +
                                                   "freq= %s, " +
                                                   "use_in_fit= %s, " +
                                                   "counterjet= %s, " +
                                                   "rating= %s" + " WHERE source='" + source + "' AND observer='" +
                                                   observer + "' AND method='UV' AND stokes='I' AND epoch='" +
                                                   row['date'] +"' AND id='" + str(row["component_number"]) + "'")

            cursor.execute(update_query, (source, row["date"], row['component_number'], row["flux"],row["radius"], row["theta"],
                                          row["size"],row["ratio"],"0","I",observer,"UV",freq,"1","-1","1"))
        # Commit the transaction
        mydb.commit()

def download_kinematic_from_MOJAVE(source,band,observer,password,difmap_path):

    # connect to database:
    mydb = sql.connect(
        host="mojavedb.mpifr-bonn.mpg.de",
        user="agn",
        password=password,
        database="galaxies"
    )

    comp_table = "components"
    bandname_uvf = "u"
    if band == 'K':
        comp_table = "Kcomponents"
        bandname_uvf = "k"
    elif band == 'Q':
        comp_table = "Qcomponents"
        bandname_uvf = "q"

    sql_query = "select * from {0} where source='{1}' and observer='{2}' and stokes='I' and method='UV' order by epoch,id".format(
        comp_table, source, observer)

    component_table = pd.read_sql(sql_query,con=mydb)

    component_table.to_csv("~/Desktop/PhD_Project/test.csv")

    epochs=np.unique(component_table["epoch"])

    #create folder
    os.makedirs("tmp_data",exist_ok=True)

    #extract data epoch by epoch
    for index,epoch in enumerate(epochs):
        #print .mod file based on MOJAVE data

        with open("tmp_data/modfile.mod", 'w') as file:

            for ind, row in component_table[component_table["epoch"]==epoch].iterrows():
                print("{:.6f}".format(float(row["flux"])),
                      "{:.6f}".format(float(row["dist"])),
                      "{:.3f}".format(float(row["pa"])),
                      "{:.6f}".format(float(row["size"])),
                      "{:.5f}".format(float(row["ratio"])),
                      "{:.3f}".format(float(row["cpa"])),
                      "1",
                      "{:.5e}".format(float(row["freq"]*1e9)),
                      "0",
                      file=file)

        epoch=str(epoch)
        #download .uvf file from MOJAVE database
        filename=source+"."+bandname_uvf+"."+epoch.replace("-","_")
        url = ("http://www.cv.nrao.edu/2cmVLBA/data/"+source+"/"+epoch.replace("-","_")+"/"+filename+".uvf")
        local_filename="tmp_data/"+filename+".uvf"  # Local file name to save the download

        """
        # Send a GET request to the URL
        response = requests.get(url, stream=True)

        # Check if the request was successful
        if response.status_code == 200:
            # Open a local file with write-binary mode
            with open(local_filename, 'wb') as file:
                # Write the response content to the file in chunks
                for chunk in response.iter_content(chunk_size=1024*1024):
                    file.write(chunk)
            print(f"File downloaded successfully as {local_filename}")
        else:
            print(f"Failed to download file. Status code: {response.status_code}")
        """

        os.chdir("tmp_data")
        #Initialize DIFMAP to create modelfit.fits file
        child = pexpect.spawn(difmap_path, encoding='utf-8', echo=False)
        child.expect_exact("0>", None, 2)

        def send_difmap_command(command, prompt="0>"):
            child.sendline(command)
            child.expect_exact(prompt, None, 2)

        send_difmap_command("obs " + filename+".uvf")
        send_difmap_command("select i")
        send_difmap_command("rmod modfile.mod")
        send_difmap_command("maps 1024")
        send_difmap_command("save "+filename)
        os.system("rm -rf difmap.log*")

        #now create component info dataframe
        try:
            if index==0:
                df_comp = getComponentInfo(filename+".fits")
            else:
                df_comp = pd.concat([df_comp, getComponentInfo(filename+".fits")],ignore_index=True)
        except:
            pass
        os.chdir("..")

    #now match the df_comp with the MOJAVE model associations:
    # Round the columns to 4 decimal places
    component_table['flux'] = component_table['flux'].round(3)
    component_table['dist'] = component_table['dist'].round(3)
    component_table['pa'] = component_table['pa'].round(3)

    df_comp['Flux_round'] = df_comp['Flux'].round(3)
    df_comp['radius_round'] = df_comp['radius'].round(3)
    df_comp['theta_round'] = df_comp['theta'].round(3)

    # Merge the DataFrames on the specified columns
    df_merged = pd.merge(df_comp, component_table, left_on=['Flux_round', 'radius_round', 'theta_round'], right_on=['flux', 'dist', 'pa'])

    # Add df1["id"] as df2["component_number"]
    df_merged['component_number'] = df_merged['id']

    #rename some stuff
    df_merged["flux"] = df_merged["Flux"]
    df_merged["maj"] = df_merged["Major_axis"]
    df_merged["min"] = df_merged["Minor_axis"]
    df_merged["pos"] = df_merged["PA"]
    df_merged["x"] = df_merged["Delta_x"]
    df_merged["y"] = df_merged["Delta_y"]
    df_merged['redshift'] = 0

    # Create the "is_core" column
    df_merged['is_core'] = df_merged['component_number'] == 0

    # Drop the now redundant columns from df_merged
    df_merged = df_merged.drop(columns=['Flux', 'dist', 'pa', 'id',"Flux_round","radius_round","theta_round",
                                        "Major_axis","Minor_axis","PA"])

    os.chdir("tmp_data")
    df_merged.to_csv("component_info.csv")
    os.makedirs("modelfit_files",exist_ok=True)
    os.system("mv *.fits modelfit_files")

    #TODO need to create empty clean_fits folder and empty kinematic_fit.csv and then can import the stuff with import_kinematics
    #TODO then need to implement buttons to do it.


#download_kinematic_from_MOJAVE("0506+056","Ku","MLL","VLBA2cm","/usr/local/difmap/uvf_difmap_2.5g/difmap")

