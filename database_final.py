import os
import os.path as op
import numpy as np
import pymysql
import shutil
import csv
import pandas as pd
import glob
from traits.api import HasTraits, Str, Enum, Range, Directory
from traitsui.api import View, Item, Handler
import pdb

def __escape_name(s):
    """ Code copied from internet source: formats string data from
    target_val_pair for error-free uploading into mySQL database.

    'Escape name to avoid SQL injection and keyword clashes.
    Doubles embedded backticks, surrounds the whole in backticks.
    Note: not security hardened, caveat emptor.''

    """
    return '`{}`'.format(s.replace('`', '``'))

def createimplanttable(cursor,db):
    '''Create the implant_db table. This should NOT be used except during development.
    May be called after the deltable function. '''
    # ---------------------------- create table for implant/region info ------------
    #cursor = db.cursor()
    cursor.execute( "CREATE TABLE implants ( \
                    implant_id INTEGER NOT NULL AUTO_INCREMENT, \
                    animal_id VARCHAR(255) NOT NULL, \
                    experiment_id VARCHAR(255), \
                    species VARCHAR(255), \
                    strain VARCHAR(255), \
                    genotype VARCHAR(255), \
                    sex VARCHAR(255), \
                    animal_dob VARCHAR(20), \
                    region VARCHAR(255) NOT NULL, \
                    daqsys VARCHAR(255), \
                    nchan SMALLINT NOT NULL, \
                    probenum TINYINT NOT NULL, \
                    chan_range VARCHAR(255) NOT NULL, \
                    n_implant_sites TINYINT, \
                    implant_date VARCHAR(255), \
                    expt_start VARCHAR(255), \
                    expt_end VARCHAR(255), \
                    surgeon VARCHAR(10), \
                    implant_coordinates VARCHAR(255) NOT NULL, \
                    electrode VARCHAR(255), \
                    headstage VARCHAR(255),\
                    PRIMARY KEY(implant_id) ) "     )

    #cursor.execute("ALTER TABLE implants ADD COLUMN implant_id INTEGER NOT NULL AUTO_INCREMENT PRIMARY KEY FIRST")
    print('Created table "implant_db" in the {} database'.format(db.db))

def createrestarttable(cursor, db):
    cursor.execute( "CREATE TABLE restarts ( \
                restart_id MEDIUMINT UNSIGNED PRIMARY KEY AUTO_INCREMENT,\
                implant_id INTEGER,\
                start_day VARCHAR(255), \
                end_day VARCHAR(255),\
                save_loc VARCHAR(255), \
                FOREIGN KEY(implant_id) REFERENCES implants(implant_id)\
                )")

def createclusterstable(cursor,db):
    '''Create the clusters table. This should NOT be used except during development.
    May be called after the deltable function. '''

    cursor.execute( "CREATE TABLE clusters ( \
                    cluster_id MEDIUMINT UNSIGNED PRIMARY KEY AUTO_INCREMENT,\
                    implant_id INTEGER, \
                    restart_id MEDIUMINT UNSIGNED,\
                    quality TINYINT, \
                    neg_pos_t DOUBLE, \
                    half_width DOUBLE, \
                    slope_falling DOUBLE, \
                    mean_amplitude DOUBLE, \
                    fr DOUBLE, \
                    cluster_idx INT, \
                    duration DOUBLE, \
                    clustering_t0 VARCHAR(255), \
                    tracklinks VARCHAR(255), \
                    folder_location VARCHAR(255),\
                    FOREIGN KEY(implant_id) REFERENCES implants(implant_id),\
                    FOREIGN KEY(restart_id) REFERENCES restarts(restart_id)\
                    )" )

    # add a column for cluster barcodes and make it the primary key and make it first.
    #cursor.execute("ALTER TABLE clusters ADD COLUMN barcode DOUBLE NOT NULL AUTO_INCREMENT PRIMARY KEY FIRST")

    # make foreign keys in the clusters table
    #cursor.execute("ALTER TABLE clusters ADD FOREIGN KEY(implant_id) REFERENCES implant_db(implant_id) ")
    print('Created table "clusters" in the {} database'.format(db.db))

def __implantgui():
    ''' This function asks a user to input the implant/region info about each electrode used in chronic electrophys recordings. The function will write a csv file into the folder that contains the relevant dataset. The csv file will be uploaded into the lab's mysql database into the implant_db table in the clusters database.'''

    # set up dictionaries for GUI field enumeration
    specieslist = {

        'Unknown'   : ['unknown'],
        'mouse'     : ['m'],
        'rat'       : ['r'],
    }

    surgeonlist = {

        'EAB'   : ['lit'],
        'KBH'   : ['kbh'],
        'LIT'   : ['lit'],
        'SCF'   : ['scf'],
        'SJB'   : ['sjb'],
        'CAF'   : ['caf']
    }

    sexlist = {
        'unknown' : ['unknown'],
        'female': ['f'],
        'male'  : ['m'],

    }

    strainlist = {

        'Unknown'       : ['unknown'],
        'Long Evans'    : ['le'],
        'Sprague Dawley': ['sd'],
        'C57 Black 6'   : ['c57b6'],

    }

    genotypelist = {

        'Unknown'   : ['unknown'],
        'WT'        : ['wt'],
        'Dnmt3a'    : ['dnmt3a'],
        'nf1'       : ['nf1'],
        'P301S/E4'  : ['p301se4'],
        'PV-cre'    : ['pvcre'],
        'Shank3'    : ['shank3'],
        'VIP-cre'   : ['vipcre'],

    }

    regionlist = {

        'Unknown'           : ['unknown'],
        'Basal Forebrain'   : ['basalforebrain'],
        'BLA'               : ['bla'],
        'CA1'               : ['ca1'],
        'CeA'               : ['cea'],
        'DG'                : ['dg'],
        'Endo ctx'          : ['endoctx'],
        'Ento ctx'          : ['entoctx'],
        'HC'                : ['hc'],
        'LGN'               : ['lgn'],
        'M1'                : ['m1'],
        'm2'                : ['m2'],
        'OFC'               : ['OFC'],
        'perirh ctx'        : ['perirh ctx'],
        'NAc'               : ['nac'],
        'RSC'               : ['rsc'],
        'SCN'               : ['scn'],
        'S1'                : ['s1'],
        'S2'                : ['s2'],
        'Subiculum'         : ['subiculum'],
        'V1m'               : ['v1m'],
        'V1b'               : ['v1b'],
        'V2'                : ['v2'],

    }


    daqsyslist = {

        'Unknown'   : ['unknown'],
        'ecube'     : ['ecube'],
        'intan'     : ['intan']

    }

    tflist = {

        'yes' : [1],
        'no'  : [0]

    }

    electrodelist = {

        'Unknown'       : ['unknown'],
        'tetrode'       : ['tetrode'],
        'stereotrode'   : ['stereotrode'],
        'single wire'   : ['single'],
        'carbon'        : ['carbon'],
        'MIT silicon'   : ['mitsilicon'],
        'UCLA silicon'  : ['uclasilicon'],
        'neuropixel'    : ['neuropixel'],
    }

    headstagelist = {

        'Unknown' : ['unknown'],
        'intan16' : ['intan16'],
        'intan32' : ['intan32'],
        'intan64' : ['intan64'],
        'HS64'    : ['hs64'],
        'HS640'   : ['hs640'],

    }

    class implantupload(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
        """
        masterpath      = Directory(os.getcwd())
        surgeon         = Str ('initials, ex. ABC')
        animalid        = Str ('ex. ABC12345')
        experiment_id   = Str ('ex. 0010101')
        species         = Enum(list(specieslist.keys())[0], list(specieslist.keys()))
        strain          = Enum(list(strainlist.keys())[0],list(strainlist.keys()))
        genotype        = Enum(list(genotypelist.keys())[0],list(genotypelist.keys()))
        sex             = Enum(list(sexlist.keys())[0],list(sexlist.keys()))
        animal_dob      = Str ('MMDDYYYY')
        region          = Enum(list(regionlist.keys())[0],list(regionlist.keys()))
        implantcoord    = Str('ex. -2.3,0.4')
        daqsys          = Enum(list(daqsyslist.keys())[0],list(daqsyslist.keys()))
        nchan           = Range(low = 1, high = 640)
        probenum        = Range(low = 1, high = 10)
        chanrange       = Str ('ex. 65-128')
        nsites          = Range(low = 1, high = 20)
        implant_date    = Str ('MMDDYYYY')
        exptstart       = Str ('MMDDYYYY')
        exptend         = Str ('MMDDYYYY')
        electrode_type  = Enum(list(electrodelist.keys())[0],list(electrodelist.keys()))
        hstype          = Enum(list(headstagelist.keys())[0],list(headstagelist.keys()))

        view = View(
            Item(name='masterpath',label='Directory'),
            Item(name = 'surgeon'),
            Item(name = 'animalid'),
            Item(name = 'experiment_id'),
            Item(name = 'species'),
            Item(name = 'strain'),
            Item(name = 'genotype'),
            Item(name = 'sex'),
            Item(name = 'animal_dob'),
            Item(name = 'region' ),
            Item(name = 'implantcoord'),
            Item(name = 'daqsys'),
            Item(name = 'electrode_type'),
            Item(name = 'hstype'),
            Item(name = 'nchan'),
            Item(name = 'probenum'),
            Item(name = 'chanrange'),
            Item(name = 'nsites'),
            Item(name = 'implant_date'),
            Item(name = 'exptstart'),
            Item(name = 'exptend'),
            title = 'Implant Information.',
            buttons = ['OK'],
            resizable = True,
            scrollable = True,

        )

    # Create the GUI:
    igui = implantupload()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        igui.configure_traits()

    return igui

def submit_implant(g,cursor,db):
    '''SUBMIT_IMPLANT Takes the data structure output from the GUI __implantgui
        and writes the data contents into the clusteringdb in the implant_db
        table. This will automatically call the clustercrawl function to crawl
        through the directory selected by the user (via the GUI), extract the
        cluster information, and write those data into the clusters table.

        Inputs:
            G: this is the output structure from __implantgui.
            CURSOR: Database cursor.
            DB: Database connection.

        Outputs:
            uniqueid,
            g.changroup,
            g.masterpath
            '''
    # convert the binary fields to 0 and 1 integers for proper formatting
    # d = np.array([g.videobin, g.lightbin, g.soundbin,g.swbin])
    # d = [int(i) for i in d == 'yes']

    # create a dictionary of each of the targets (names) and the correspoding data
    target_val_pair = {
            "animal_id": g.animalid.upper(),
            "experiment_id" : g.experiment_id,
            "species" : g.species,
            "sex" : g.sex,
            "animal_dob" : g.animal_dob,
            "region" : g.region,
            "strain" : g.strain,
            "genotype" : g.genotype,
            "daqsys" : g.daqsys,
            "nchan" : g.nchan,
            "probenum"  : g.probenum,
            "chan_range" : g.chanrange,
            "n_implant_sites" : g.nsites,
            "implant_date" : g.implant_date,
            "expt_start" : g.exptstart,
            "expt_end" : g.exptend,
            "surgeon" : g.surgeon,
            "implant_coordinates" : g.implantcoord,
            "electrode" : g.electrode_type,
            "headstage" :g.hstype
            }

    # convert dictionary target and value information into tuples
    targets = tuple( [*target_val_pair] )
    #values  = tuple( [*target_val_pair.values()] )
    # automatically rewrite the target names into the proper format for mysql
    cols = ', '.join(map(__escape_name, targets))  # assumes the keys are *valid column names*.
    placeholders = ', '.join(['%({})s'.format(name) for name in targets])
    # submit to the implants_db table.
    print(g.nchan)
    query = 'INSERT INTO implants ({}) VALUES ({})'.format(cols, placeholders)

    cursor.execute(query, target_val_pair)
    uniqueid = cursor.lastrowid

    db.commit()
    print('Added implant information to the implant_db table in the clusteringdb database.')

    # add the folder location and the implant barcode ID (unique, generated on
    # commit) to the dictionary
    target_val_pair.update({'location':g.masterpath, 'implant_id':uniqueid})

    # write to a pandas dataframe, use this to write to a .csv file easily.
    df = pd.DataFrame.from_dict(data=target_val_pair, orient='index')
    fn = g.masterpath + '/' + g.animalid + '_' + g.region + '_' + str(g.nsites) + '_sites.csv'
    pd.DataFrame.from_dict(data=target_val_pair, orient='index').to_csv(fn, header=False)
    print('Wrote implant information to .csv file  {}'.format(fn))

    return uniqueid, g.probenum, g.masterpath
    # This information should be sent to the clustercrawl function. clustercrawl
    # will automatically calculate/detect cluster metadata by implant (channel
    # group) and block (time) and write to another table in the database.