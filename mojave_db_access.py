import mysql.connector as sql
import pandas as pd
import numpy as np
import requests
import os
import pexpect
from graph_generator import getComponentInfo, get_date
import subprocess
from datetime import datetime, timedelta
import math

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

        #test if the date matches the date from the database (can be different by +/-1 days)
        def get_mojave_date(date,source,iter=0):
            response = requests.get("https://www.cv.nrao.edu/2cmVLBA/data/"+source+"/"+date.replace("-","_")+"/", stream=True)
            # Check if the request was successful
            if response.status_code == 200:
                return date
            elif (response.status_code == 404) and (iter <= 4):
                # in this case, probably the MOJAVE date is out of sync +/-1 day with the .uvf date, let's try to match it....
                iter += 1
                delta_t = (-1) ** iter * timedelta(days=iter)
                # try a different date, they might be misaligned by +/-1 day
                date_obj = datetime.strptime(date, "%Y-%m-%d")
                day_new = date_obj + delta_t
                # reset epoch to new date
                date = day_new.strftime("%Y-%m-%d")
                return get_mojave_date(date, source, iter=iter)
            else:
                return date
        row["date"]=get_mojave_date(row["date"],source)

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

def download_kinematic_from_MOJAVE(source,band,observer,password,difmap_path,foldername="tmp_data"):

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
    component_table["epoch"]=component_table["epoch"].astype("str")
    epochs=np.unique(component_table["epoch"])

    #create folder
    os.system("rm -rf "+foldername)
    os.makedirs(foldername,exist_ok=True)

    #extract data epoch by epoch
    for index,epoch in enumerate(epochs):

        #download .uvf file from MOJAVE database
        filename=source+"."+bandname_uvf+"."+epoch.replace("-","_")
        url = ("https://www.cv.nrao.edu/2cmVLBA/data/"+source+"/"+epoch.replace("-","_")+"/"+filename+".uvf")
        local_filename=foldername+"/"+filename+".uvf"  # Local file name to save the download


        #this function downloads MOJAVE .uvf files for a given epoch, if it is not findable it tries searching in +/-2days
        #since sometimes there are misalignments between the epoch and the date in the .uvf file
        def download_file(url,local_filename,filename,epoch_in,iter=0):
            response=requests.get(url, stream=True)
            # Check if the request was successful
            if response.status_code == 200:
                # Open a local file with write-binary mode
                with open(local_filename, 'wb') as file:
                    # Write the response content to the file in chunks
                    for chunk in response.iter_content(chunk_size=1024*1024):
                        file.write(chunk)
                print(f"File downloaded successfully as {local_filename}")
                return filename
            elif (response.status_code == 404) and (iter <=4):
                #in this case, probably the MOJAVE date is out of sync +/-1 day with the .uvf date, let's try to match it....
                iter+=1
                delta_t=(-1)**iter*timedelta(days=iter)
                #try a different date, they might be misaligned by +/-1 day
                date_obj=datetime.strptime(epoch_in,"%Y-%m-%d")

                day_new=date_obj+delta_t
                #reset epoch to new date
                epoch_in = day_new.strftime("%Y-%m-%d")
                #update the filenames
                filename = source + "." + bandname_uvf + "." + epoch_in.replace("-", "_")
                url = ("https://www.cv.nrao.edu/2cmVLBA/data/" + source + "/" + epoch_in.replace("-",
                                                                                              "_") + "/" + filename + ".uvf")
                local_filename = foldername + "/" + filename + ".uvf"  # Local file name to save the download
                return download_file(url,local_filename,filename,epoch_in,iter=iter)
            else:
                print(f"Failed to download file. Status code: {response.status_code}")
                return filename

        filename=download_file(url,local_filename,filename,epoch)
        try:
            epoch_new=get_date(foldername+"/"+filename+".uvf")
            epochs[index]=epoch_new
            component_table.loc[component_table["epoch"] == epoch,"epoch"] = epoch_new
            epoch=epoch_new
        except:
            pass


        # print .mod file based on MOJAVE database data
        with open(foldername + "/modfile.mod", 'w') as file:

            for ind, row in component_table[component_table["epoch"] == epoch].iterrows():
                print("{:.6f}".format(float(row["flux"])),
                      "{:.6f}".format(float(row["dist"])),
                      "{:.3f}".format(float(row["pa"])),
                      "{:.6f}".format(float(row["size"])),
                      "{:.5f}".format(float(row["ratio"])),
                      "{:.3f}".format(float(row["cpa"])),
                      "1",
                      "{:.5e}".format(float(row["freq"] * 1e9)),
                      "0",
                      file=file)

        os.chdir(foldername)
        #Initialize DIFMAP to create modelfit.fits file
        child = pexpect.spawn(difmap_path+"/difmap", encoding='utf-8', echo=False, cwd=foldername)
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

        #Here we need to modify the epoch date again to match the .fits file (for some reason .uvf and .fits date can be different)
        try:
            epoch_new = get_date(foldername + "/" + filename + ".fits")
            epochs[index] = epoch_new
            component_table.loc[component_table["epoch"] == epoch, "epoch"] = epoch_new
            epoch = epoch_new
        except:
            pass

        #now create component info dataframe
        try:
            df_comp = pd.concat([df_comp, getComponentInfo(foldername+"/"+filename+".fits")],ignore_index=True)
        except:
            try:
                df_comp = getComponentInfo(foldername+"/"+filename+".fits")
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
    df_merged = pd.merge(df_comp, component_table, left_on=["Date",'Flux_round', 'radius_round', 'theta_round'], right_on=['epoch','flux', 'dist', 'pa'])

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

    os.chdir(foldername)
    df_merged.to_csv("component_info.csv")
    os.makedirs("modelfit_files",exist_ok=True)
    os.system("mv *.fits modelfit_files")
    os.makedirs("uvf_files", exist_ok=True)
    os.system("mv *.uvf uvf_files")
    os.makedirs("clean_fits",exist_ok=True)
    os.system("cp component_info.csv kinematic_fit.csv")
    output = subprocess.check_output("pwd",shell=True)
    output = output.decode("utf-8").strip()
    os.chdir("..")

def query_models(source,password):

    # connect to database:
    mydb = sql.connect(
        host="mojavedb.mpifr-bonn.mpg.de",
        user="agn",
        password=password,
        database="galaxies"
    )

    bands=["Ku","K","Q"]
    options=[]
    for i,table in enumerate(["components","Kcomponents","Qcomponents"]):
        sql_query="select observer from {0} where source='{1}' group by observer".format(table,source)
        df=pd.read_sql(sql_query,con=mydb)
        for ind,obs in enumerate(df["observer"]):
            options.append(obs+" "+bands[i])

    return options

def check_password(password):
    # try connecting to database:
    try:
        mydb = sql.connect(
            host="mojavedb.mpifr-bonn.mpg.de",
            user="agn",
            password=password,
            database="galaxies"
        )
    except:
        return False
    return True

#download_kinematic_from_MOJAVE("0506+056","Ku","MLL","VLBA2cm","/usr/local/difmap/uvf_difmap_2.5g/difmap")

