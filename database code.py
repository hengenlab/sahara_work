# to make a new database
import pymysql
db = pymysql.connect(user = 'root', password = 'N3tw0rks#!')

cursor = db.cursor()
cursor.execute('CREATE DATABASE clusteringdb')

# if it says you didnt select a database
db.select_db(db) #idk why this works, seems redundant but it worked
cursor = db.cursor() # to get the new cursor


# directories to test shit
# /Volumes/HlabShare/mbttest/inputmod/ti
# /Volumes/HlabShare/Clustering_Data/EAB00040/ts/
# /Volumes/HlabShare/Clustering_Data/SCF00008/52919-10lam/t
# Documents/GitHub/LabTestFiles/new_clustering
# /Volumes/HlabShare/Clustering_Data/EAB00048/Clustered_EAB48_06-10_06-11/final/

query = "DESCRIBE clusters"
query = 'SELECT * FROM clusters'
cursor.execute(query)
records = cursor.fetchall()

### TO SET UP DB!! ONLY FOR DEV
import pymysql
db = pymysql.connect(user = 'root', password = 'N3tw0rks#!')

cursor = db.cursor()
cursor.execute('CREATE DATABASE clusteringdb')
cursor, db = connectclusterdb('root', 'N3tw0rks#!')
createimplanttable(cursor, db)
createclusterstable(cursor, db)



### CODE TO GIVE EVERYONE
import musclebeachtools as mbt
cursor, db = mbt.database.connectclusterdb('root', 'N3tw0rks#!') # this connects you to the database, might not need this but good to do it anyway
mbt.database.upload_implant('root','N3tw0rks#!') # make sure the path is right or youre in the right directory





    # query = "  SELECT implant_db.animal_id FROM clusters JOIN implant_db ON clusters.implant_id = implant_db.implant_id WHERE clusters.barcode IN ({}) ".format(res2)
    #
    # query = 'INSERT INTO implant_db ({}) VALUES ({})'.format(cols, placeholders)





# # "where" search:
    # query = "SELECT barcode FROM clusters WHERE mean_amplitude > 10 AND clusters.implant_id = implant_db.implant_id )  "
    # records = cursor.fetchall()
    #
    #
    #
    #
    # (SELECT salesman_id
    #  FROM salesman
    #  WHERE name='Paul Adam');
    #
    #
    #
    # # SHOW  RECORDS IN TABLE - - - - - - - - - - -
    # retreive one column only
query = "SELECT quality FROM clusters"
# retreive all columns
query = "SELECT * FROM implant_db"
#retrieve some columns
query = "SELECT quality, mean_amplitude FROM clusters"
# show the column np_samples
query = "DESCRIBE clusters"
## getting records from the table
cursor.execute(query)
## fetching all records from the 'cursor' object
records = cursor.fetchall()
## Showing the data
for record in records:
    print(record)
    # - - - - - - - - - - - - - - - - - - - - - - - - -


