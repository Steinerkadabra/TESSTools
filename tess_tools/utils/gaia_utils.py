from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import sys
import numpy as np
import matplotlib.pyplot as plt
from astroquery.mast import Catalogs

from astroquery.gaia import Gaia

def add_gaia_figure_elements(tpf, magnitude_limit=18,targ_mag=10.):
    """Make the Gaia Figure Elements"""
    # Get the positions of the Gaia sources
    c1 = SkyCoord(tpf.ra, tpf.dec, frame='icrs', unit='deg')
    # Use pixel scale for query size
    pix_scale = 4.0  # arcseconds / pixel for Kepler, default
    if tpf.mission == 'TESS':
        pix_scale = 21.0
    # We are querying with a diameter as the radius, overfilling by 2x.
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    result = Vizier.query_region(c1, catalog=["I/345/gaia2"],
                                 radius=Angle(np.max(tpf.shape[1:]) * pix_scale, "arcsec"))
    no_targets_found_message = ValueError('Either no sources were found in the query region '
                                          'or Vizier is unavailable')
    too_few_found_message = ValueError('No sources found brighter than {:0.1f}'.format(magnitude_limit))
    if result is None:
        raise no_targets_found_message
    elif len(result) == 0:
        raise too_few_found_message
    result = result["I/345/gaia2"].to_pandas()
    result = result[result.Gmag < magnitude_limit]
    if len(result) == 0:
        raise no_targets_found_message
    radecs = np.vstack([result['RA_ICRS'], result['DE_ICRS']]).T
    coords = tpf.wcs.all_world2pix(radecs, 0.5) ## TODO, is origin supposed to be zero or one?
    year = ((tpf.time[0].jd - 2457206.375) * u.day).to(u.year)
    pmra = ((np.nan_to_num(np.asarray(result.pmRA)) * u.milliarcsecond/u.year) * year).to(u.arcsec).value
    pmdec = ((np.nan_to_num(np.asarray(result.pmDE)) * u.milliarcsecond/u.year) * year).to(u.arcsec).value
    result.RA_ICRS += pmra
    result.DE_ICRS += pmdec

    # Gently size the points by their Gaia magnitude
    sizes = 128.0 / 2**(result['Gmag']/targ_mag)#64.0 / 2**(result['Gmag']/5.0)
    one_over_parallax = 1.0 / (result['Plx']/1000.)
    r = (coords[:, 0]+tpf.column,coords[:, 1]+tpf.row,result['Gmag'])

    return r,result

# Plot orientation
def plot_orientation(tpf):
	"""
    Plot the orientation arrows

    Returns
    -------
    tpf read from lightkurve

	"""
	mean_tpf = np.mean(tpf.flux,axis=0)
	nx,ny = np.shape(mean_tpf)
	x0,y0 = tpf.column+int(0.8*nx),tpf.row+int(0.2*nx)
	# East
	tmp =  tpf.get_coordinates()
	ra00, dec00 = tmp[0][0][0][0], tmp[1][0][0][0]
	ra10,dec10 = tmp[0][0][0][-1], tmp[1][0][0][-1]
	theta = np.arctan((dec10-dec00)/(ra10-ra00))
	if (ra10-ra00) < 0.0: theta += np.pi
	#theta = -22.*np.pi/180.
	x1, y1 = 1.*np.cos(theta), 1.*np.sin(theta)
	plt.arrow(x0,y0,x1,y1,head_width=0.2,color='white')
	plt.text(x0+1.5*x1,y0+1.5*y1,'E',color='white')
	# North
	theta = theta +90.*np.pi/180.
	x1, y1 = 1.*np.cos(theta), 1.*np.sin(theta)
	plt.arrow(x0,y0,x1,y1,head_width=0.2,color='white')
	plt.text(x0+1.5*x1,y0+1.5*y1,'N',color='white')



def get_gaia_data(ra, dec):
    """
    Get Gaia parameters

    Returns
    -------
    RA, DEC
    """
    # Get the positions of the Gaia sources
    c1 = SkyCoord(ra, dec, frame='icrs', unit='deg')
    # We are querying with a diameter as the radius, overfilling by 2x.
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    result = Vizier.query_region(c1, catalog=["I/345/gaia2"],
                                 radius=Angle(10., "arcsec"))
    try:
        result = result["I/345/gaia2"]
    except:
        print('Not in Gaia DR2. If you know the Gaia ID and Gmag, try the options --gid and --gmag.')
        print('Exiting without finishing...')
        sys.exit()

    no_targets_found_message = ValueError('Either no sources were found in the query region '
                                          'or Vizier is unavailable')
    too_few_found_message = ValueError('No sources found closer than 1 arcsec to TPF coordinates')
    if result is None:
        raise no_targets_found_message
    elif len(result) == 0:
        raise too_few_found_message

    if len(result)>1:
        dist = np.sqrt((result['RA_ICRS']-ra)**2 + (result['DE_ICRS']-dec)**2)
        idx = np.where(dist == np.min(dist))[0][0]
        return result[idx]['Source'], result[idx]['Gmag']
    else:
        return result[0]['Source'], result[0]['Gmag']

def get_gaia_data_from_tic(tic):
    '''
    Get Gaia parameters

    Returns
    -----------------------
    GaiaID, Gaia_mag
    '''
    # Get the Gaia sources
    result = Catalogs.query_object('TIC'+tic, radius=.005, catalog="TIC")
    IDs = result['ID'].data.data
    k = np.where(IDs == tic)[0][0]
    GAIAs = result['GAIA'].data.data
    Gaiamags = result['GAIAmag'].data.data

    GAIA_k = GAIAs[k]
    Gaiamag_k = Gaiamags[k]

    if GAIA_k == '':
        GAIA_k = np.nan
    return GAIA_k, Gaiamag_k


def get_coord(tic):
	"""
	Get TIC corrdinates

	Returns
	"""
	try:
		catalog_data = Catalogs.query_object(objectname="TIC"+tic, catalog="TIC")
		ra = catalog_data[0]["ra"]
		dec = catalog_data[0]["dec"]
		return ra, dec
	except:
		print("ERROR: No gaia ID found for this TIC")
          

def dr3_from_dr2(dr2ID):
    query_dr3fromdr2 = "select dr3_source_id from gaiadr3.dr2_neighbourhood where dr2_source_id = "+dr2ID
    job = Gaia.launch_job(query=query_dr3fromdr2)
    dr3_ids = job.results['dr3_source_id'].value.data
    if len(dr3_ids) == 1:
        myid = dr3_ids[0]
    else:
        print("\t WARNING! There are more than one DR3 ids for this DR2 ID, assuming the first one...")
        myid = dr3_ids[0]

    return myid

def get_gaia_data_from_simbad(dr2ID):
    # simb = Simbad.query_object('Gaia DR2 '+dr2ID)
    # simbid = Simbad.query_objectids('Gaia DR2 '+dr2ID)
    # if simbid == None:
    #     print("ERROR: TIC not found in Simbad as Gaia DR2 "+str(dr2ID))
    # ids = np.array(simbid['ID'].data).astype(str)
    # myid = [id for id in ids if 'DR3' in id]
    # if len(myid) == 0:
    #     myid = [id for id in ids if 'DR2' in id]
    # myid = myid[0].split(' ')[2]

    myid = dr3_from_dr2(dr2ID)
    query2 = "SELECT \
             TOP 1 \
             source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag\
             FROM gaiadr3.gaia_source\
             WHERE source_id = "+str(myid)+" \
             "
    job = Gaia.launch_job_async(query2)
    gmag = job.get_results()['phot_g_mean_mag'].data[0]

    return myid,gmag
