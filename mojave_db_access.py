import mysql.connector as sql
import pandas as pd

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