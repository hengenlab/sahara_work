import os
import numpy as np
import pymysql
import pandas as pd
import base64
from pandas.io.sql import DatabaseError
from traits.api import HasTraits, Str, Enum, Range, Directory, Bool, File, Int, List, Date, Float
from traitsui.api import View, Item, Handler, Action, VGroup, FileEditor, CheckListEditor, RangeEditor
from traitsui.menu import OKCancelButtons

def __escape_name(s):
    """ Code copied from internet source: formats string data from
    target_val_pair for error-free uploading into mySQL database.

    'Escape name to avoid SQL injection and keyword clashes.
    Doubles embedded backticks, surrounds the whole in backticks.
    Note: not security hardened, caveat emptor.''

    """
    return '`{}`'.format(s.replace('`', '``'))


def __get_conditionals(g):
    """
        Returns a string that can be used in a query to pull out
        information that fits these conditions
    """
    a = []
    s = []
    m = []
    z = []
    if g.animal != 'any':
        a = [__format_conditional('animals.animal_name', g.animal)]
    if g.site != 'any':
        s = [__format_conditional('probes.region', g.site)]
    if g.manipulation != 'any':
        m = [__format_conditional('restarts.manipulation', g.manipulation)]
    if g.genotypes != 'any':
        z = [__format_conditional('animals.genotype', g.genotypes)]

    all_conds = a + s + m + z
    return ' AND '.join(all_conds)


def __get_joins(list, cond):
    '''
    Returns a string that can be used in a query
    the 'join' string since the queries can only join
    the tables its using
    '''
    check = ' '.join(list) + ' ' + cond
    str = ''
    if 'probes' in check:
        str = str + 'JOIN probes ON animals.animal_id = probes.animal_id '
    if 'restarts' in check:
        str = str + 'JOIN restarts ON animals.animal_id = restarts.animal_id '
    if 'clusters' in check:
        str = str + 'JOIN clusters ON animals.animal_id = clusters.animal_id '
    return str


def format_query(table, s):
    '''
    joins the table and the field together
    '''
    return f'{table}.{s}'


def __format_conditional(table, con):
    '''
    joins the table + field and the condition it must meet
    '''
    return f'{table} = "{con}"'


def __check_existance(name, cursor, db):
    '''
    returns true if that animal exists in the animals table
    '''
    q = f'SELECT EXISTS(SELECT * FROM animals WHERE animal_name = "{name}")'
    cursor.execute(q)
    result = np.asarray(cursor.fetchall()).flatten()[0]
    return result

def __check_existance_restart(name, start, cursor):
    q = f'SELECT EXISTS(SELECT * FROM restarts WHERE animal_id = {name} AND start_day = "{start}")'
    cursor.execute(q)
    result = np.asarray(cursor.fetchall()).flatten()[0]
    return result

def connectclusterdb(user, pwd):
    ''' CONNECTCLUSTERDB. Connect to the clusteringdb database.
    Inputs:
        USER: username
        PWD: password

    Outputs:
        DB.CURSOR: "A database cursor is an identifier associated with a group
            of rows. It is, in a sense, a pointer to the current row in a buffer.
            You must use a cursor in the following cases: Statements that return
            more than one row of data from the database server: A SELECT statement
            requires a select cursor."
        DB: Database connection.
     '''
    # connect to the clustering database

    db = pymysql.connect(user   = user,
                                passwd  = pwd,
                                host    = "localhost",
                                database= "clusteringdb")


    return db.cursor(), db


# #create table functions not up to date, see MYSQL workbench query files
# def __createanimaltable(cursor, db):
#     cursor.execute("CREATE TABLE animals ( \
#                     animal_id INTEGER NOT NULL AUTO_INCREMENT, \
#                     animal_name VARCHAR(255) NOT NULL, \
#                     species VARCHAR(255), \
#                     strain VARCHAR(255), \
#                     genotype VARCHAR(255), \
#                     sex VARCHAR(255), \
#                     animal_dob VARCHAR(20), \
#                     num_chan SMALLINT NOT NULL, \
#                     num_sites TINYINT NOT NULL, \
#                     implant_date VARCHAR(255), \
#                     surgeon VARCHAR(10), \
#                     electrode VARCHAR(255), \
#                     headstage VARCHAR(255),\
#                     PRIMARY KEY(animal_id) ) ")
#
#
# def __createprobetable(cursor,db):
#     '''Create the implant_db table. This should NOT be used except during development.
#     May be called after the deltable function. '''
#     # ---------------------------- create table for implant/region info ------------
#     #cursor = db.cursor()
#     cursor.execute( "CREATE TABLE probes ( \
#                     probe_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
#                     animal_id INTEGER NOT NULL, \
#                     probe_num SMALLINT NOT NULL, \
#                     region VARCHAR(255), \
#                     chan_range VARCHAR(255), \
#                     FOREIGN KEY(animal_id) REFERENCES animals(animal_id)")
#
#     #cursor.execute("ALTER TABLE implants ADD COLUMN implant_id INTEGER NOT NULL AUTO_INCREMENT PRIMARY KEY FIRST")
#     print('Created table "implant_db" in the {} database'.format(db.db))
#
#
# def __createrestarttable(cursor, db):
#     cursor.execute( "CREATE TABLE restarts ( \
#                 restart_id MEDIUMINT UNSIGNED PRIMARY KEY AUTO_INCREMENT,\
#                 animal_id INTEGER,\
#                 start_day VARCHAR(255), \
#                 end_day VARCHAR(255),\
#                 save_loc VARCHAR(255), \
#                 manipulations VARCHAR(255), \
#                 FOREIGN KEY(animal_id) REFERENCES animals(animal_id)\
#                 )")
#
#
# def __createclusterstable(cursor,db):
#     """Create the clusters table. This should NOT be used except during development.
#     May be called after the deltable function. """
#
#     cursor.execute( "CREATE TABLE clusters ( \
#                     cluster_id MEDIUMINT UNSIGNED PRIMARY KEY AUTO_INCREMENT,\
#                     animal_id INTEGER, \
#                     probe_id INTEGER, \
#                     restart_id MEDIUMINT UNSIGNED,\
#                     quality TINYINT, \
#                     neg_pos_t DOUBLE, \
#                     half_width DOUBLE, \
#                     slope_falling DOUBLE, \
#                     mean_amplitude DOUBLE, \
#                     fr DOUBLE, \
#                     cluster_idx INT, \
#                     duration DOUBLE, \
#                     clustering_t0 VARCHAR(255), \
#                     tracklinks VARCHAR(255), \
#                     folder_location VARCHAR(255),\
#                     FOREIGN KEY(animal_id) REFERENCES animals(animal_id),\
#                     FOREIGN KEY(probe_id) REFERENCES probes(probe_id), \
#                     FOREIGN KEY(restart_id) REFERENCES restarts(restart_id)\
#                     )" )
#
#     # add a column for cluster barcodes and make it the primary key and make it first.
#     #cursor.execute("ALTER TABLE clusters ADD COLUMN barcode DOUBLE NOT NULL AUTO_INCREMENT PRIMARY KEY FIRST")
#
#     # make foreign keys in the clusters table
#     #cursor.execute("ALTER TABLE clusters ADD FOREIGN KEY(implant_id) REFERENCES implant_db(implant_id) ")
#     print('Created table "clusters" in the {} database'.format(db.db))


def __probegui(animal_name, title):
    """
    GUI class for uploading to the probe database
    return the gui object with all the information from the
    user

    parameters:
    animal_name - current animal you're uploading
    title - probe number you're on

    return:
    gui object
    """
    regionlist = ['Unknown', 'Basal Forebrain', 'BLA', \
        'CA1','CeA','DG','Endo ctx','Ento ctx','HC',\
        'LGN','M1','m2','OFC','perirh ctx','NAc', \
        'RSC','SCN','S1','S2','Subiculum','V1m',\
        'V1b','V2']


    class probeupload(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
        """
        animal          = Str (animal_name[0])
        probe_num       = Int (title)
        chanrange       = Str ('ex. 65-128')
        region          = Enum(regionlist)
        ap              = Float()
        ml              = Float()
        dv              = Float()
        exit            = Bool()

        view = View(

            Item(name = 'animal', style = 'readonly'),
            Item(name = 'probe_num', style = 'readonly'),
            Item(name = 'chanrange'),
            Item(name = 'region'),
            Item(name = 'ap', label = 'AP'),
            Item(name = 'ml', label = 'ML'),
            Item(name = 'dv', label = 'DV'),
            Item(name = 'exit', label = 'EXIT'),

            title = f'Probe {title} Information.',
            buttons = ['OK'],
            resizable = True,
            scrollable = True,
        )

    # Create the GUI:
    pgui = probeupload()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        pgui.configure_traits()

    return pgui


def __animalgui():
    """
    GUI for uploading to the animal table

    Very top level gui, doesn't take any other information into it

    returns: gui object with information from user
    """

    # set up dictionaries for GUI field enumeration
    specieslist = ['Unknown', 'Mouse', 'Rat']

    surgeonlist = ['CAF', 'EAB', 'KBH', 'SCF', 'SJB']

    sexlist = ['Unknown', 'Female', 'Male']

    strainlist = ['Unknown', 'Long Evans', 'Sprague Dawley', 'C57 Black 6']

    genotypelist = ['Unknown', 'WT', 'Dnmt3a', 'nf1', 'P301s/E4', 'E4', 'PVcre', 'Shank3', 'AI', 'VIPcre']

    regionlist = ['Unknown', 'Basal Forebrain', 'BLA', \
                  'CA1', 'CeA', 'DG', 'Endo ctx', 'Ento ctx', 'HC', \
                  'LGN', 'M1', 'm2', 'OFC', 'perirh ctx', 'NAc', \
                  'RSC', 'SCN', 'S1', 'S2', 'Subiculum', 'V1m', \
                  'V1b', 'V2']

    daqsyslist = ['Unknown', 'ecube', 'intan']

    electrodelist = ['Unknown', 'tetrode', 'sterotrode', 'single wire', 'carbon', 'MIT silicon', 'UCLA silicon', 'neuropixel']

    headstagelist = ['Unknown', 'intan16', 'intan32', 'intan64', 'HS64', 'HS640']

    class animalupload(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
        """
        exit = Bool()
        surgeon         = Enum(list(surgeonlist))
        animalid        = Str ('ex. ABC12345')
        species         = Enum(list(specieslist))
        strain          = Enum(list(strainlist))
        genotype        = Enum(list(genotypelist))
        sex             = Enum(list(sexlist))
        animal_dob      = Date()
        daqsys          = Enum(list(daqsyslist))
        nchan           = Range(low = 1, high = 640)
        numprobes        = Range(low = 1, high = 10)
        nregions         = Range(low = 1, high = 20)
        implant_date    = Date()
        electrode_type  = Enum(list(electrodelist))
        hstype          = Enum(list(headstagelist))
        active          = Bool()
        sac_date        = Date()

        view = View(
            Item(name = 'exit', label = 'EXIT'),
            Item(name = 'surgeon'),
            Item(name = 'animalid'),
            Item(name = 'implant_date'),


            Item(name = 'species'),
            Item(name = 'strain'),
            Item(name = 'genotype'),
            Item(name = 'sex'),
            Item(name = 'animal_dob'),


            Item(name = 'daqsys'),
            Item(name = 'electrode_type'),
            Item(name = 'hstype'),


            Item(name = 'nchan'),
            Item(name = 'numprobes'),
            Item(name = 'nregions'),

            Item(name = 'active'),
            Item(name = 'sac_date'),


            title = 'Implant Information.',
            buttons = ['OK'],
            resizable = True,
            scrollable = True,
        )

    # Create the GUI:
    agui = animalupload()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        agui.configure_traits()

    return agui


def __restartgui(current_animals, ogstart_day, ogend_day, ogcameras):
    manipulationslist = ['none', 'MD', 'Food Dep', 'Crickets', 'Cocaine', 'Visual Stim', 'Auditory Dep', 'Operant Task', 'Sprinkles', 'OPTO', 'DREADDs', 'Lick Port']

    class restartupload(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
        """
        masterpath      = Directory(os.getcwd())
        animalid        = Enum(list(current_animals))
        start_day       = Date(ogstart_day)
        end_day         = Date(ogend_day)
        cameras         = Str(ogcameras)
        manipulations   = Enum(list(manipulationslist))
        exit = Bool()

        view = View(
            Item(name = 'exit', label = 'EXIT'),
            Item(name = 'masterpath', label = 'Directory'),
            Item(name = 'animalid'),
            Item(name = 'start_day'),
            Item(name = 'end_day'),
            Item(name = 'cameras'),
            Item(name = 'manipulations'),
            title = 'Restart Information.',
            buttons = ['OK'],
            resizable = True,
            scrollable = True,
        )

    # Create the GUI:
    rgui = restartupload()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        rgui.configure_traits()

    return rgui


def __clustergui(current_animals, current_restarts):
    """
    Gui for uploading to the clustering table

    parameters:
    current animals - current animals in the database
    current restarts - current restarts in the database

    returns:
    cluster gui object with all the information filled out by the user
    """
    class clusterupload(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
        """
        masterpath      = File(editor = FileEditor())
        animal          = Enum(list(current_animals))
        restart         = Enum(list(current_restarts))
        probe_num       = Enum(list(np.arange(1, 9)))
        time_frame      = Str('ex. 0-12')
        exit = Bool()
        view = View(
            Item(name = 'exit', label = 'EXIT'),
            Item(name = 'masterpath', label="Clustering File"),
            Item(name = 'animal', label = 'Animal'),
            Item(name = 'restart', label = 'Restart'),
            Item(name = 'probe_num', label = 'Probe Number'),
            Item(name = 'time_frame', label = 'Time Frame Clustered'),
            title = 'Cluster Information.',
            buttons = ['OK'],
            resizable = True,
            scrollable = True,
        )

    # Create the GUI:
    cgui = clusterupload()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        cgui.configure_traits()

    return cgui


def __sac_gui(current_animals):
    """
    Small gui to fill out saccing information about currently running animals

    parameters:
    current_animals - list of currently running animals ("alive" = True)

    returns:
    gui object with info the user filled out
    """
    class sac(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
        """
        animal        = Enum(list(current_animals))
        sac_date      = Date()
        exit = Bool()
        view = View(
            Item(name = 'exit', label = 'EXIT'),
            Item(name = 'animal'),
            Item(name = 'sac_date'),
            title = 'Sac Information.',
            buttons = ['OK'],
            resizable = True,
            scrollable = True,
        )

    # Create the GUI:
    sgui = sac()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        sgui.configure_traits()

    return sgui

def __dlcsws_gui(current_animals, current_restarts):


    class dlcsws(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
        """
        animal          = Enum (list(current_animals))
        restart         = Enum (list(current_restarts))
        sleep_wake_extracted = Bool()
        sleep_wake_scored = Bool()
        DLC_motion_location = Directory()
        DLC_features = Str()
        exit = Bool()
        view = View(

            Item(name = 'animal'),
            Item(name = 'restart'),
            "_",
            Item(name = 'sleep_wake_extracted', label = 'Sleep Wake Extraction'),
            Item(name = 'sleep_wake_scored', label = 'Sleep Wake Scored'),
            "_",
            Item(name = 'DLC_motion_location', label = 'DLC Motion Location'),
            Item(name = 'DLC_features', label = 'DLC Training Features'),
            "_",
            Item(name = 'exit', label = 'EXIT'),
            title = f'DLC + SWS Editor',
            buttons = ['OK'],
            resizable = True,
            scrollable = True,
        )

    # Create the GUI:
    pgui = dlcsws()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        pgui.configure_traits()

    return pgui


def submit_clusters(g, cursor, db):
    """
    This is the function that submits the clustering outputs to the database itself

    goes through the path you gave it and runs brief stats on all the clusters then
    submits each cluster to the database

    returns: the last id from the entered clusters
    """
    samplerate=25000
    #get implant and restart info
    q = f'SELECT animal_id FROM clusteringdb.animals WHERE animal_name = "{g.animal}"'
    cursor.execute(q)
    animal_id = np.asarray(cursor.fetchall()).flatten()[0]

    q2 = f'SELECT restart_id FROM clusteringdb.restarts WHERE animal_id = {animal_id} AND start_day = "{g.restart}"'
    q3 = f'SELECT probe_id FROM clusteringdb.probes WHERE animal_id = {animal_id} AND probe_num = {g.probe_num}'
    try:
        cursor.execute(q2)
        restart_id = np.asarray(cursor.fetchall()).flatten()[0]

        cursor.execute(q3)
        probe_id = np.asarray(cursor.fetchall()).flatten()[0]

    except IndexError:
        print('There is no restart or probe connected to that animal in the database. Please add it before you submit a clustering job.')
        return None

    #first gotta load the cells
    try:
        cells = np.load(g.masterpath, allow_pickle = True)
        count = 0
        for cell in cells:
            neg_pos_t = (cell.waveforms[cell.waveforms.argmin():-1].argmax()) * 1/samplerate
            half = np.sum( cell.waveforms < 0.5 * cell.waveforms.min() ) * (1/samplerate)
            falling = (cell.waveforms[28] - cell.waveforms[24]) / (4 * (1 / samplerate) * 1000)
            fr = len(cell.spike_time)/cell.end_time

            clust_stats = {
                'animal_id': int(animal_id),
                'restart_id' : int(restart_id),
                'probe_id' : int(probe_id),
                'quality' : int(cell.quality),
                'neg_pos_t' : float(neg_pos_t),
                'half_width' : float(half),
                'slope_falling' : float(falling),
                'mean_amplitude' : int(cell.mean_amplitude),
                'fr' : float(fr),
                'cluster_idx' : int(cell.clust_idx),
                'duration' : float(cell.start_time - cell.end_time),
                'clustering_t0' : g.restart,
                'tracklinks' : 'lol',
                'folder_location' : g.masterpath,
                'time_frame' : g.time_frame
                }

            targets = tuple([*clust_stats])
            cols = ', '.join(map(__escape_name, targets))
            placeholders = ', '.join(['%({})s'.format(name) for name in targets])
            query = 'INSERT INTO clusters ({}) VALUES ({})'.format(cols, placeholders)
            cursor.execute(query, clust_stats)
            uniqueid = cursor.lastrowid
            count += 1
    except ValueError as err:
        print("There was an error uploading your clusters. Make sure your clustering file is clean.")
        print(f"ERROR: {str(err)}")
        return
    db.commit()
    print(f'Submitted {count} clusters to the clusters table. cluster_ids {uniqueid-count} - {uniqueid}')
    return uniqueid


def submit_restart(g, cursor, db):
    """
    This is the function that actually submits the restart information from the
    gui object to the restarts table

    it links from the animal to the aniamls table and then submits all the
    info from the user

    returns the unique id of the restart entry
    """
    query1 = f'SELECT animal_id FROM clusteringdb.animals WHERE animal_name = "{g.animalid}"'
    cursor.execute(query1)
    animal_id = np.asarray(cursor.fetchall()).flatten()[0]

    target_val_pair = {
        "animal_id": int(animal_id),
        "start_day": g.start_day,
        "end_day": g.end_day,
        "cameras" : g.cameras,
        "save_loc":g.masterpath,
        "manipulations":g.manipulations
    }
    # convert dictionary target and value information into tuples
    targets = tuple( [*target_val_pair] )

    exists = __check_existance_restart(animal_id, g.start_day, cursor)
    if exists:
        print("This restart already exists for this animal, if you want to edit it go write an edit query.")
        return None
    else:
        #values  = tuple( [*target_val_pair.values()] )
        # automatically rewrite the target names into the proper format for mysql
        cols = ', '.join(map(__escape_name, targets))  # assumes the keys are *valid column names*.
        placeholders = ', '.join(['%({})s'.format(name) for name in targets])
        query = 'INSERT INTO restarts ({}) VALUES ({})'.format(cols, placeholders)
        cursor.execute(query, target_val_pair)
        uniqueid = cursor.lastrowid

        print(f'Submitted {g.start_day} restart with id {uniqueid} to restarts table')
        db.commit()
        return uniqueid


def submit_probe(g, cursor, db):
    '''This function takes the info from the gui
    and submits it to the probe table of the database

    links to the animal table by animal id and submits probe information from user

    returns: unique id of the probe entry
    '''

    q = f'SELECT animal_id FROM clusteringdb.animals WHERE animal_name = "{g.animal}"'
    cursor.execute(q)
    animal_id = np.asarray(cursor.fetchall()).flatten()[0]

    # convert the binary fields to 0 and 1 integers for proper formatting
    # d = np.array([g.videobin, g.lightbin, g.soundbin,g.swbin])
    # d = [int(i) for i in d == 'yes']

    # create a dictionary of each of the targets (names) and the correspoding data


    target_val_pair = {
        "animal_id": int(animal_id),
        "probe_num": int(g.probe_num),
        "region": g.region,
        "ap": g.ap,
        "ml":g.ml,
        "dv":g.dv,
        "chan_range": g.chanrange
    }

    # convert dictionary target and value information into tuples
    targets = tuple([*target_val_pair])
    # values  = tuple( [*target_val_pair.values()] )
    # automatically rewrite the target names into the proper format for mysql
    cols = ', '.join(map(__escape_name, targets))  # assumes the keys are *valid column names*.
    placeholders = ', '.join(['%({})s'.format(name) for name in targets])
    query = 'INSERT INTO probes ({}) VALUES ({})'.format(cols, placeholders)

    cursor.execute(query, target_val_pair)
    uniqueid = cursor.lastrowid

    db.commit()
    print('Added probe information to the probes table in the clusteringdb database.')

    # add the folder location and the implant barcode ID (unique, generated on
    # commit) to the dictionary
    # target_val_pair.update({'location':g.masterpath, 'implant_id':uniqueid})
    #
    # # write to a pandas dataframe, use this to write to a .csv file easily.
    # df = pd.DataFrame.from_dict(data=target_val_pair, orient='index')
    # fn = g.masterpath + '/' + g.animalid + '_' + g.region + '_' + str(g.nsites) + '_sites.csv'
    # pd.DataFrame.from_dict(data=target_val_pair, orient='index').to_csv(fn, header=False)
    # print('Wrote implant information to .csv file  {}'.format(fn))

    return uniqueid


def submit_animal(g,cursor,db):
    '''SUBMIT_ANIMAL Takes the data structure output from the GUI __animalgui
        and writes the data contents into the lab database in the animals
        table.

        Inputs:
            G: this is the output structure from __implantgui.
            CURSOR: Database cursor.
            DB: Database connection.

        Outputs:
            uniqueid of the animal entry
            number of probes from this animal
            '''
    # convert the binary fields to 0 and 1 integers for proper formatting
    # d = np.array([g.videobin, g.lightbin, g.soundbin,g.swbin])
    # d = [int(i) for i in d == 'yes']
    # create a dictionary of each of the targets (names) and the correspoding data
    target_val_pair = {
            "animal_name": g.animalid.upper(),
            "species" : g.species,
            "sex" : g.sex,
            "animal_dob" : str(g.animal_dob),
            "strain" : g.strain,
            "genotype" : g.genotype,
            "daqsys" : g.daqsys,
            "num_chan" : g.nchan,
            "num_probes"  : g.numprobes,
            "num_regions" : g.nregions,
            "implant_date" : str(g.implant_date),
            "alive" : g.active,
            "sac_date" : str(g.sac_date),
            "surgeon" : g.surgeon,
            "electrode" : g.electrode_type,
            "headstage" :g.hstype
            }

    # convert dictionary target and value information into tuples
    exists = __check_existance(g.animalid.upper(), cursor, db)
    if exists:
        print('That animal already exists, if you want to edit it please write an edit query, see the help docs for the search funciton.')
        return None, None
    else:
        if g.active:
            target_val_pair.pop('sac_date')
        targets = tuple( [*target_val_pair] )
        #values  = tuple( [*target_val_pair.values()] )
        # automatically rewrite the target names into the proper format for mysql
        cols = ', '.join(map(__escape_name, targets))  # assumes the keys are *valid column names*.
        placeholders = ', '.join(['%({})s'.format(name) for name in targets])
        # submit to the implants_db table.

        query = 'INSERT INTO animals ({}) VALUES ({})'.format(cols, placeholders)

        cursor.execute(query, target_val_pair)
        uniqueid = cursor.lastrowid

        db.commit()
        print('Added implant information to the aniamls table in the lab_bd database.')
        print(f'submitted a {g.nchan} channel animal')
        # add the folder location and the implant barcode ID (unique, generated on
        # commit) to the dictionary
        # target_val_pair.update({'location':g.masterpath, 'implant_id':uniqueid})
        #
        # # write to a pandas dataframe, use this to write to a .csv file easily.
        # df = pd.DataFrame.from_dict(data=target_val_pair, orient='index')
        # fn = g.masterpath + '/' + g.animalid + '_' + g.region + '_' + str(g.nsites) + '_sites.csv'
        # pd.DataFrame.from_dict(data=target_val_pair, orient='index').to_csv(fn, header=False)
        # print('Wrote implant information to .csv file  {}'.format(fn))

        return uniqueid, int(g.numprobes)


    # This information should be sent to the clustercrawl function. clustercrawl
    # will automatically calculate/detect cluster metadata by implant (channel
    # group) and block (time) and write to another table in the database.


def submit_sac(g, cursor, db):
    q = f'SELECT animal_id FROM clusteringdb.animals WHERE animal_name = "{g.animal}"'
    cursor.execute(q)
    animal_id = np.asarray(cursor.fetchall()).flatten()[0]

    query = f'UPDATE animals SET alive = 0, sac_date = "{str(g.sac_date)}" WHERE animal_id = {animal_id}'

    cursor.execute(query)

    db.commit()
    print('Updated the sac date in the animals table.')

    return

def submit_dlcsws(g, cursor, db):
    q = f'SELECT animal_id FROM clusteringdb.animals WHERE animal_name = "{g.animal}"'
    cursor.execute(q)
    animal_id = np.asarray(cursor.fetchall()).flatten()[0]

    q = f'SELECT restart_id FROM clusteringdb.restarts WHERE animal_id = {animal_id} AND start_day = "{g.restart}"'
    try:
        cursor.execute(q)
        restart_id = np.asarray(cursor.fetchall()).flatten()[0]
    except:
        print("There are no restarts with that animal and start day in the database")
        return

    if g.DLC_features == '' or g.DLC_motion_location == '':
        query = f'UPDATE restarts SET sleep_wake_extracted = {g.sleep_wake_extracted}, sleep_wake_scored = {g.sleep_wake_scored} WHERE restart_id = {restart_id}'
    else:
        query = f'UPDATE restarts SET sleep_wake_extracted = {g.sleep_wake_extracted}, sleep_wake_scored = {g.sleep_wake_scored}, DLC_motion_location = "{g.DLC_motion_location}", DLC_features = "{g.DLC_features}" WHERE restart_id = {restart_id}'
    cursor.execute(query)
    db.commit()
    print('Updated the restart info')



def upload_animal(user,pwd):
    '''Top level function that calls the gui and submitting function for
    the animals table

    returns: unique animal id and number of probes, to be used in the main looped function
    '''
    #connect to the clusteringdb
    cursor, db = connectclusterdb (user, pwd)
    # call the implant info GUI and collect relevant information.
    g = __animalgui()
    if g.exit:
        print('--EXITING--')
        return 0,0
    # format and pass information from GUI to the implant_db table
    uniqueid, numprobes = submit_animal(g,cursor,db)
    return uniqueid, numprobes


def upload_probe(user, pwd, animal_id, probenum):
    """
    Top level probe uploading function that calls the gui and the
    submitting function for a probe into the probes table
    """
    cursor, db = connectclusterdb(user, pwd)

    query = f'SELECT animal_name FROM clusteringdb.animals WHERE animal_id = "{animal_id}"'
    cursor.execute(query)
    animal_name = np.asarray(cursor.fetchall()).flatten()

    g = __probegui(animal_name, probenum)
    if g.exit:
        print('--EXITING--')
        return
    uniqueid = submit_probe(g, cursor, db)


def upload_restart(user, pwd):
    """
    top level function for the restart table, calls the gui and the submitting function
    """
    cursor, db = connectclusterdb(user, pwd)

    query = 'SELECT animal_name FROM clusteringdb.animals'
    cursor.execute(query)
    current_implants = np.asarray(cursor.fetchall()).flatten()

    query = 'SELECT * FROM clusteringdb.restarts ORDER BY restart_id DESC LIMIT 1'
    result = pd.read_sql(query, db)

    if result.empty:
        sday = None
        eday = None
        cs = None
    else:
        sday = result.start_day[0]
        eday = result.end_day[0]
        cs = result.cameras[0]

    g = __restartgui(current_implants, sday, eday, cs)

    if g.exit:
        print('--EXITING--')
        return
    uniqueid = submit_restart(g, cursor, db)


def upload_dlc_sws_edit(user, pwd):
    cursor, db = connectclusterdb(user, pwd)

    query = 'SELECT animal_name FROM clusteringdb.animals'
    cursor.execute(query)
    current_animals = np.asarray(cursor.fetchall()).flatten()

    query = 'SELECT start_day FROM clusteringdb.restarts'
    cursor.execute(query)
    current_restarts = np.unique(np.asarray(cursor.fetchall()).flatten())

    g = __dlcsws_gui(current_animals, current_restarts)
    if g.exit:
        print('--EXITING--')
        return
    submit_dlcsws(g, cursor, db)


def upload_clusters(user, pwd):
    """
    top level function for submitting clustering outputs to the database
    """
    cursor, db = connectclusterdb(user, pwd)

    query = 'SELECT animal_name FROM clusteringdb.animals'
    cursor.execute(query)
    current_animals = np.asarray(cursor.fetchall()).flatten()

    query = 'SELECT start_day FROM clusteringdb.restarts'
    cursor.execute(query)
    current_restarts = np.unique(np.asarray(cursor.fetchall()).flatten())

    g = __clustergui(current_animals, current_restarts)
    if g.exit:
        print('--EXITING--')
        return
    uniqueid = submit_clusters(g, cursor, db)


def sac_animal(user, pwd):
    """top level function for submitting an animals sac data"""
    cursor, db = connectclusterdb(user, pwd)
    query = 'SELECT animal_name FROM clusteringdb.animals WHERE alive = 1'
    cursor.execute(query)
    current_animals = np.asarray(cursor.fetchall()).flatten()

    g = __sac_gui(current_animals)
    if g.exit:
        print('--EXITING--')
        return
    submit_sac(g, cursor, db)


def __top_gui():
    """
    Gui for the top most function
    this is where the user will decide what action they want to take with the database
    """

    class myHandler(Handler):
        def _close_me(self, info):
            info.ui.dispose()


    class toplevel(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
        """
        CLOSE = False

        action = Enum(["DONE",'submit a surgery', 'submit a restart', 'submit clusters', 'sac', 'DLC/SWS', 'QUERY'])
        close_button = Action(name="Close", label = "close", action= '_close_me')
        view = View(
            Item(name = 'action'),
            handler = myHandler(),
            title = 'Database access',
            buttons = ['OK'],
            resizable = True,
            scrollable = True
        )

    # Create the GUI:
    top_gui = toplevel()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        top_gui.configure_traits()

    return top_gui


def use_the_database():
    """
    This is the function that will be called from the command line to interact with the database
    consists of a loop that continues going until a user specifies that they're done using the database

    returns:
    None - if there is no dataframe requested
    df - a pandas dataframe of your query outputs as long as your query is sucessful




    QUERY CRASH COURSE:
           There are 4 tables in this database. animals, probes, restarts, and clusters. If you want to pull information
           from more than 1 table then you need to JOIN them. Every single table is
           joined on the .animal_id field, as that is the linking quality for all of our
           data.

           To get information from the database, you SELECT what field you want from
           what tables. Then you JOIN the tables based ON a certain field. Then you
           decide WHERE the database should select your data.

           If you're unsure what fields are available for each table, take a peak at the PNG stored in the
           database folder of HlabShare. It's a visual representation of the database with all the fields and connections between the tables.

    QUERY SAMPLES:
        simple:
            SELECT * FROM animals
            --- pulls all the entries from the animal table

            SELECT * FROM animals WHERE surgeon = 'CAF'
            --- pulls all fields from any entry where Clayton was the surgeon.

            SELECT animal_name, animal_genotype FROM animals
            ---  Only outputs the animals name and genotype for every entry in the animals table.

        moderate:
            SELECT animals.animal_name, animals.genotype FROM animals
            JOIN probes ON animals.animal_id = probes.animal_id
            WHERE probes.region = 'CA1'
            --- Outputs the animals name and genotype for animals that had a probe in CA1

            SELECT animals.animal_name, restarts.start_day, restarts.save_loc FROM animals
            JOIN restarts ON restarts.animal_id = animals.animals_id
            WHERE restarts.manipulation = 'MD'
            --- Outputs the animal name, restart start day and save location for any restart + animal
            that had MD during that block

        complex:
            SELECT animals.animal_name, clusters.cluster_idx, clusters.folder_location FROM animals
            JOIN clusters ON clusters.animal_id = animals.animal_id
            JOIN probes ON probes.animal_id = animals.animal_id
            WHERE probes.region = 'M1' AND clusters.quality = 1
            --- Outputs the animal name, cluster idx and clustering output location of any cluster
            that was from M1 and had a quality of 1.
            --- This way you can have a list of the idxs and the location of the clustering files to load

            SELECT animals.animal_name, probes.region, animals.genotype FROM animals
            JOIN probes ON animals.animal_id = probes.animal_id
            WHERE probes.region = 'CA1' OR probes.region = 'M1'
            --- pulls animal name, probe region, and the genotype for any animal that had a
            probe in CA1 or M1

            SELECT animals.animal_name, probes.region, probes.probe_num, animals.genotype, restarts.start_day, restarts.save_loc FROM animals
            JOIN probes ON animals.animal_id = probes.animal_id
            JOIN restarts ON animals.animal_id = restarts.animal_id
            WHERE restarts.sleep_wake_scored = 1
            --- pulls animal name, probe region, probe number, animal genotype, restart, and raw data save location
            for all restarts that have been sleep scored

        most complex: (AND USEFUL)
            SELECT DISTINCT animals.animal_name, restarts.start_day, clusters.folder_location, 'YES' as clustered FROM animals
            JOIN clusters ON animals.animal_id = clusters.animal_id
            JOIN restarts ON clusters.restart_id = restarts.restart_id
            --- returns information for all restarts that have been clustered, and adds a column that says they were clustered, so the output makes sense.
            also only returns 1 row if there are repeats (there will be hundreds since its pulling from the clusters table). Notice the way
            the joins are structured to make sure its not attaching restarts to animals that have not been clustered.

            SELECT DISTINCT animals.animal_name, restarts.start_day, restarts.save_loc, restarts.restart_id, 'NO' as clustered FROM animals
            JOIN restarts ON animals.animal_id = restarts.animal_id
            WHERE restarts.restart_id NOT IN (SELECT clusters.restart_id FROM clusters)
            --- returns information about any restart that HASNT been clustered yet. Essentially pulls from all restarts that exist in the
            restart table but don't exist in the clustering table.



    You probably won't need to get any more complex than that since this is mainly just an organizational
    system. But if you want to get fancy with how the output is sorted, grouped, or do any kind of calculation
    in the query - Google is your friend.

    """
    t = []
    with open('/media/HlabShare/db_creds.txt') as f:
        for line in f:
            t.append(line.rstrip())

    pwd = base64.b64decode(t[0]).decode()
    user = base64.b64decode(t[1]).decode()

    alive = True
    to_return = None
    while alive:
        g = __top_gui()
        if g.action == 'DONE':
            print("thanks for visiting")
            alive = False

        elif g.action == 'submit a surgery':
            animal_id, numprobes = upload_animal(user, pwd)
            if animal_id is not None:
                for p in range(numprobes):
                    upload_probe(user, pwd, animal_id, p+1)
        elif g.action == 'submit a restart':
            upload_restart(user, pwd)
        elif g.action == 'submit clusters':
            upload_clusters(user, pwd)
        elif g.action == 'sac':
            sac_animal(user, pwd)
        elif g.action == 'QUERY':
            to_return = SEARCH(user, pwd)
        elif g.action == 'DLC/SWS':
            upload_dlc_sws_edit(user, pwd)


    return to_return


def SEARCH(user, pwd):
    """
    Searching function

    gets the outputs from the different search guis and prepares and executes queries based off
    the users input

    returns:
    df - a pandas dataframe with the results of the specified queries
    none - if the query wasn't sucessful it returns none

    """
    cursor, db = connectclusterdb(user, pwd)

    query = 'SELECT animal_name FROM clusteringdb.animals'
    cursor.execute(query)
    current_animals = np.asarray(cursor.fetchall()).flatten()
    if len('current_animals')==0:
        current_animals = ['any']
    else:
        current_animals = np.insert(current_animals, 0, 'any')

    query = 'SELECT manipulations FROM clusteringdb.restarts'
    cursor.execute(query)
    current_manipulations = np.unique(np.asarray(cursor.fetchall()).flatten())
    if len(current_manipulations) == 0:
        current_manipulations = ['any']
    else:
        current_manipulations= np.insert(current_manipulations, 0, 'any')


    query = 'SELECT region FROM clusteringdb.probes'
    cursor.execute(query)
    current_sites = np.unique(np.asarray(cursor.fetchall()).flatten())
    if len(current_sites)==0:
        current_sites = ['any']
    else:
        current_sites = np.insert(current_sites, 0, 'any')

    query = 'SELECT genotype FROM clusteringdb.animals'
    cursor.execute(query)
    current_genos = np.unique(np.asarray(cursor.fetchall()).flatten())
    if len(current_genos) == 0:
        current_genos = ['any']
    else:
        current_genos = np.insert(current_genos, 0, 'any')

    g = __search_gui(current_sites, current_manipulations, current_genos, current_animals, 'main_view')
    if g.exit:
        print('--EXITING--')
        return None
    df = None
    if g.own_query != ' ':
        q = g.own_query

        if 'ALTER' in q.upper() or 'UPDATE' in q.upper():
            cont = input("You're about to update something in the database, are you sure?  y/n ")
            if cont != 'y':
                return None
        if 'DROP' in q.upper():
            print('Not a chance. If you think you need to drop a table go ask Sahara.')
            return None
        try:
            df = pd.read_sql(q, db)
        except DatabaseError:
            print("Your query didn't work ----- sorry ")
            print(DatabaseError)

        if df.empty:
            print('Sorry. There is nothing matching that query in the database')
        else:
            print("QUERY SUCCESSFUL -- dataframe will be returned when you exit")
            print(df)

    else:
        conditionals = __get_conditionals(g)
        if len(conditionals) == 0:
            print('You have to give it a parameter. If you want all the data we have, go look at the simple queries')
            return df
        if g.descriptions:
            d = __desc_gui()
            total_map = d.animal_descs + d.restart_descs + d.probe_descs + d.cluster_descs

            cols = ', '.join(total_map)
            joins = __get_joins(total_map, conditionals)

            q = f'SELECT DISTINCT {cols} ' \
                f'FROM animals ' \
                f'{joins}' \
                f'WHERE {conditionals}'

            df = pd.read_sql(q, db)
            if df.empty:
                print("Sorry, there is nothing matching that description in the database")
            else:
                print(df)

        if g.clusters:
            if 'clusters' not in joins:
                joins += "JOIN clusters ON animals.animal_id = clusters.animal_id "
            q = 'SELECT animals.animal_name, clusters.cluster_idx, clusters.folder_location ' \
                'FROM animals ' \
                f'{joins}' \
                f'WHERE {conditionals}'
            df = pd.read_sql(q, db)
            if df.empty:
                print('Sorry, there are no clusters matching that description')
            else:
                print(df)
    return df


def __desc_gui():
    """
    gui for the descriptives output

    consists of many checkboxes where you can decide what output you want to see
    from the database if you chose that you want descriptives
    """
    class descgui(HasTraits):
        animal_info_value_list = ['animal_name', 'surgeon', 'animal_dob', 'species', 'strain',
                                  'genotype', 'sex', 'num_chan', 'num_sites', 'num_probes',
                                  'implant_date', 'electrode', 'headstage', 'daqsys', 'alive', 'sac_date']
        animal_info_value_list = list(map(format_query, np.repeat('animals', len(animal_info_value_list)), animal_info_value_list))
        animal_info_key_list = ['Animal Name', 'Surgeon', 'Animal DOB', 'Species', 'Strain',
                                'Genotype', 'Sex','# Channels', '# Sites', '# Probes',
                                'Implant Date', 'Electrode Type', 'Headstage', 'Daqsys', 'Active?', 'Sac Date']

        animal_info = list(zip(animal_info_value_list, animal_info_key_list))

        probe_info_value_list = ['region', 'probe_num', 'chan_range', 'implantcoords']
        probe_info_value_list = list(map(format_query, np.repeat('probes', len(probe_info_value_list)), probe_info_value_list))

        probe_info_key_list = ['Region', 'Probe Number', 'Channel Range', 'Implant Coordinates']
        probe_info = list(zip(probe_info_value_list, probe_info_key_list))

        restart_info_value_list = ['start_day', 'end_day', 'save_loc', 'manipulations', 'sleep_wake_extracted', 'sleep_wake_scored', 'DLC_motion_location', 'DLC_features']
        restart_info_key_list = ['Start Date', 'End Date', 'Save Location', 'Manipulations', 'Sleep Wake Extracted', 'Sleep Wake Scored', 'DLC Motion Save Location', 'DLC Training Features']
        restart_info_value_list = list(map(format_query, np.repeat('restarts', len(restart_info_value_list)), restart_info_value_list))

        restart_info = list(zip(restart_info_value_list, restart_info_key_list))

        cluster_info_value_list = ['quality', 'neg_pos_t', 'half_width', 'slope_falling',
                                   'mean_amplitude', 'fr', 'cluster_idx', 'duration', 'folder_location']
        cluster_info_key_list = ['Quality', 'Neg-Pos Time', 'Half Width', 'Slope Falling', 'Mean Amplitude',
                                 'Mean Firing Rate', 'Cluster Idx', 'Duration', 'Save Location']
        cluster_info_value_list = list(map(format_query, np.repeat('clusters', len(cluster_info_value_list)), cluster_info_value_list))
        cluster_info = list(zip(cluster_info_value_list, cluster_info_key_list))

        animal_descs = List(animal_info, editor = CheckListEditor(values = animal_info, cols = 4))
        probe_descs = List(probe_info, editor = CheckListEditor(values = probe_info, cols = 4))
        restart_descs = List(restart_info, editor = CheckListEditor(values = restart_info, cols = 4))
        cluster_descs = List(cluster_info, editor = CheckListEditor(values = cluster_info, cols = 4))

        traits_view = View(
            Item(name = 'animal_descs', label = "Animal Descriptors", style = 'custom'),
            '_',
            Item(name = 'probe_descs', label = "Probe Descriptors", style = 'custom'),
            '_',
            Item(name = 'restart_descs', label = "Restart Descriptors", style = 'custom'),
            '_',
            Item(name = 'cluster_descs', label = "Cluster Descriptors", style = 'custom'),
            title = 'Return Descriptors',
            buttons = ['OK'],
            resizable = True,
            scrollable = True
        )

    sgui = descgui()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        sgui.configure_traits()

    return sgui


def __search_gui(current_sites, current_manipulations, current_genotypes, all_animals, view):
    """
    gui meant for the search engine

    lets the user choose what conditions they want their output to meet.

    There is also an option to write your own query.

    """
    class searchgui(HasTraits):

        site = Enum(list(current_sites))
        manipulation = Enum(list(current_manipulations))
        genotypes = Enum(list(current_genotypes))
        animal = Enum(list(all_animals))
        clusters = Bool()
        descriptions = Bool()
        own_query = Str(' ')
        exit = Bool()

        main_view = View(
            Item(name = 'site'),
            Item(name = 'manipulation'),
            Item(name = 'genotypes'),
            Item(name = 'animal'),
            Item(name = 'clusters'),
            Item(name = 'descriptions'),
            Item(name = 'own_query'),
            Item(name = 'exit', label = 'EXIT'),
            title = 'SEARCH',
            buttons = ['OK'],
            resizable = True,
            scrollable = True
        )


    # Create the GUI:
    sgui = searchgui(view = view)

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        sgui.configure_traits(view = view)

    return sgui