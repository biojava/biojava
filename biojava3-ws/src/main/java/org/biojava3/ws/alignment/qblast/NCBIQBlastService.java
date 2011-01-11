/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava3.ws.alignment.qblast;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.HashMap;

import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;

import org.biojava3.ws.alignment.RemotePairwiseAlignmentProperties;
import org.biojava3.ws.alignment.RemotePairwiseAlignmentService;
import org.biojava3.ws.alignment.RemotePairwiseAlignmentOutputProperties;

/**
 * NCBIQBlastService - A simple way of submitting BLAST request to the QBlast
 * service at NCBI.
 * 
 * <p>
 * NCBI provides a Blast server through a CGI-BIN interface. NCBIQBlastService
 * simply encapsulates an access to it by giving users access to get/set methods
 * to fix sequence, program and database as well as advanced options.
 * </p>
 * 
 * <p>
 * The philosophy behind NCBIQBlastService is to disconnect submission of
 * Blast requests from collection of Blast results. This is done so to allow a
 * user to submit multiple Blast requests while allowing recovery of the reports
 * at a later time.
 * </p>
 * 
 * <p>
 * Presently, only blastall programs are accessible. blastpgp and megablast are
 * high-priorities.
 * </p>
 * 
 * @author Sylvain Foisy, Diploide BioIT
 * @since Biojava 3
 * 
 */
public class NCBIQBlastService implements RemotePairwiseAlignmentService {

	private static String baseurl = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi";
	private URL aUrl;
	private URLConnection uConn;
	private OutputStreamWriter fromQBlast;
	private BufferedReader rd;
	
	private String email = "anonymous@biojava.org";
	private String tool = "biojavax";
	 

	private String seq = null;
	private String prog = null;
	private String db = null;
//	private String advanced = null;

	private String rid;
	private long step;
	private long start;
	private HashMap<String, Long> holder = new HashMap<String, Long>();

	/**
	 * The constructor for a QBlast service request.
	 * 
	 * @throws Exception
	 *             if the NCBI URL is unresponsive
	 * 
	 */
	public void NCBIQBlastService() throws Exception {
		try {
			aUrl = new URL(baseurl);
			uConn = setQBlastServiceProperties(aUrl.openConnection());
		}
		/*
		 * Needed but should never be thrown since the URL is static and known
		 * to exist
		 */
		catch (MalformedURLException e) {
			throw new Exception(
					"It looks like the URL for NCBI QBlast service is wrong.\n");
		}
		/*
		 * Intercept if the program can't connect to QBlast service
		 */
		catch (IOException e) {
			throw new Exception(
					"Impossible to connect to QBlast service at this time. Check your network connection.\n");
		}
	}

	/**
	 * This class is the actual worker doing all the dirty stuff related to
	 * sending the Blast request to the CGI_BIN interface. It should never be
	 * used as is but any method wanting to send a Blast request should manage
	 * to use it by feeding it the right parameters.
	 * 
	 * @param str : a <code>String</code> representation of a sequence from either of the
	 *              three wrapper methods
	 * @param rpa :a <code>RemotePairwiseAlignmentProperties</code> object
	 * @return rid : the ID of this request on the NCBI QBlast server
	 * @throws Exception if unable to connect to the NCBI QBlast service
	 */
	private String sendActualAlignementRequest(String str,
			RemotePairwiseAlignmentProperties rpa) throws Exception {

		seq = "QUERY=" + str;
		prog = "PROGRAM=" + rpa.getAlignmentOption("PROGRAM");
		db = "DATABASE=" + rpa.getAlignmentOption("DATABASE");

		if (prog == null || db == null || str == null || str.length() == 0) {
			throw new Exception(
					"Impossible to execute QBlast request. One or more of sequence|database|program has not been set correctly.\n");
		}

		String cmd = "CMD=Put&SERVICE=plain&" + seq + "&" + prog + "&"
				+ db + "&" + "FORMAT_TYPE=HTML"+"&TOOL="+getTool()+"&EMAIL="+getEmail();

		/* 
		 * This is a not so good hack to be fix by forcing key 
		 * checking in RemoteQBlastAlignmentProperties 
		 *
		 * */
		if (rpa.getAlignmentOption("OTHER_ADVANCED") != null) {
			cmd += cmd + "&" + rpa.getAlignmentOption("OTHER_ADVANCED");
		}

		try {

			uConn = setQBlastServiceProperties(aUrl.openConnection());

			fromQBlast = new OutputStreamWriter(uConn.getOutputStream());

			fromQBlast.write(cmd);
			fromQBlast.flush();

			// Get the response
			rd = new BufferedReader(new InputStreamReader(uConn
					.getInputStream()));

			String line = "";

			while ((line = rd.readLine()) != null) {
				if (line.contains("RID")) {
					String[] arr = line.split("=");
					rid = arr[1].trim();
				} else if (line.contains("RTOE")) {
					String[] arr = line.split("=");
					step = Long.parseLong(arr[1].trim()) * 1000;
					start = System.currentTimeMillis() + step;
				}
				holder.put(rid, start);
			}
		} catch (IOException e) {
			throw new Exception(
					"Can't submit sequence to BLAST server at this time.\n");
		}

		return rid;
	}

	/**
	 * This method is a wrapper that executes the Blast request via the Put
	 * command of the CGI-BIN interface with the specified parameters and a
	 * string representing the sequence. It gets the estimated time of
	 * completion by capturing the value of the RTOE variable and sets a loop
	 * that will check for completion of analysis at intervals specified by
	 * RTOE.
	 * 
	 * <p>
	 * It also capture the value for the RID variable, necessary for fetching
	 * the actual results after completion.
	 * </p>
	 * 
	 * @param str : a <code>String</code> with a sequence
	 * @param rpa : a <code>RemotePairwiseAlignmentProperties</code> object
	 * @return rid : a <code>String</code> with the request ID for this sequence
	 * @throws Exception if it is not possible to sent the BLAST command
	 */
	public String sendAlignmentRequest(String str,
			RemotePairwiseAlignmentProperties rpa) throws Exception {

		/*
		 * sending the command to execute the Blast analysis
		 */
		return rid = sendActualAlignementRequest(seq, rpa);
	}

	/**
	 * This method is a wrapper that executes the Blast request via the Put
	 * command of the CGI-BIN interface with the specified parameters and a
	 * Sequence. It gets the estimated time of completion by capturing the
	 * value of the RTOE variable and sets a loop that will check for completion
	 * of analysis at intervals specified by RTOE.
	 * 
	 * <p>
	 * It also capture the value for the RID variable, necessary for fetching
	 * the actual results after completion.
	 * </p>
	 * 
	 * @param rs :a <code>Sequence</code> object
	 * @param rpa :a <code>RemotePairwiseAlignmentProperties</code> object
	 * @return rid : a <code>String</code> with the request ID for this sequence
	 * 
	 * @throws Exception if it is not possible to sent the BLAST command
	 */
	public String sendAlignmentRequest(Sequence rs,
			RemotePairwiseAlignmentProperties rpa) throws Exception {

		seq = rs.seqString();

		return rid = sendActualAlignementRequest(seq, rpa);
	}

	/**
	 * This method is a wrapper that executes the Blast request via the Put
	 * command of the CGI-BIN interface with the specified parameters and a
	 * GenBank GID. It gets the estimated time of completion by capturing the
	 * value of the RTOE variable and sets a loop that will check for completion
	 * of analysis at intervals specified by RTOE.
	 * 
	 * <p>
	 * It also capture the value for the RID variable, necessary for fetching
	 * the actual results after completion.
	 * </p>
	 * 
	 * @param gid :an integer with a Genbank GID
	 * @param rpa :a <code>RemotePairwiseAlignmentProperties</code> object
	 * @return rid : a String with the request ID for this sequence
	 * @throws Exception if it is not possible to sent the BLAST command
	 */
	public String sendAlignmentRequest(int gid,
			RemotePairwiseAlignmentProperties rpa) throws Exception {

		seq = Integer.toString(gid);
		return rid = sendActualAlignementRequest(seq, rpa);
	}

	/**
	 * <p>
	 * This method is used only for the executeBlastSearch method to check for
	 * completion of request using the NCBI specified RTOE variable
	 * </p>
	 * 
	 * @param id : a valid request ID
	 * @param present : a representation of "now" using System.currentTimeMillis().
	 * @return a boolean value telling if the request has been completed or not.
	 */
	public boolean isReady(String id, long present) throws Exception {
		boolean isReady = false;
		String check = "CMD=Get&RID=" + id;

		if (holder.containsKey(id)) {
			/*
			 * If present time is less than the start of the search added to
			 * step obtained from NCBI, just do nothing ;-)
			 * 
			 * This is done so that we do not send zillions of requests to the
			 * server. We do the waiting internally first.
			 */
			if (present < start) {
				isReady = false;
			}
			/*
			 * If we are at least step seconds in the future from the actual
			 * call sendAlignementRequest()
			 */
			else {
				try {
					uConn = setQBlastServiceProperties(aUrl.openConnection());

					fromQBlast = new OutputStreamWriter(uConn.getOutputStream());
					fromQBlast.write(check);
					fromQBlast.flush();

					rd = new BufferedReader(new InputStreamReader(uConn
							.getInputStream()));

					String line = "";

					while ((line = rd.readLine()) != null) {
						if (line.contains("READY")) {
							isReady = true;
						} else if (line.contains("WAITING")) {
							/*
							 * Else, move start forward in time... for the next
							 * iteration
							 */
							start = present + step;
							holder.put(id, start);
						}
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		} else {
			throw new Exception("Impossible to check for request ID named "
					+ id + " because it does not exists!\n");
		}
		return isReady;
	}

	/**
	 * <p>
	 * This method extracts the actual Blast report for this request ID. It uses
	 * an object implementing the RemotePairwiseAlignmentOutputProperties
	 * interface which will specify output formatting options.
	 * </p>
	 * 
	 * @param id :a valid request ID
	 * @param rb : a <code>RemotePairwiseAlignmentOutputProperties</code> that
	 *             will specify specific output formatting commands
	 * @return an <code>InputStream</code> that can be use any way one might
	 *         desire
	 * @throws Exception if it is not possible to recover the results.
	 */
	public InputStream getAlignmentResults(String id,
			RemotePairwiseAlignmentOutputProperties rb) throws Exception {
		if (holder.containsKey(id)) {
			String srid = "CMD=Get&RID=" + id + "&"
					+ rb.getOutputOption("FORMAT_TYPE") + "&"
					+ rb.getOutputOption("ALIGNMENT_VIEW") + "&"
					+ rb.getOutputOption("DESCRIPTIONS") + "&"
					+ rb.getOutputOption("ALIGNMENTS")
					+ "&TOOL="+getTool()+"&EMAIL="+getEmail();

			try {
				uConn = setQBlastServiceProperties(aUrl.openConnection());

				fromQBlast = new OutputStreamWriter(uConn.getOutputStream());
				fromQBlast.write(srid);
				fromQBlast.flush();

				return uConn.getInputStream();

			} catch (IOException ioe) {
				throw new Exception(
						"It is not possible to fetch Blast report from NCBI at this time.\n");
			}
		} else {
			throw new Exception(
					"Impossible to get output for request ID named " + id
							+ " because it does not exists!\n");
		}
	}

	/**
	 * A simple method to check the availability of the QBlast service
	 * 
	 * @throws Exception if unable to connect to the NCBI QBlast service
	 */
	public void printRemoteBlastInfo() throws Exception {
		try {
			OutputStreamWriter out = new OutputStreamWriter(uConn
					.getOutputStream());

			out.write("CMD=Info");
			out.flush();

			// Get the response
			BufferedReader rd = new BufferedReader(new InputStreamReader(uConn
					.getInputStream()));

			String line = "";

			while ((line = rd.readLine()) != null) {
				System.out.println(line);
			}

			out.close();
			rd.close();
		} catch (IOException e) {
			throw new Exception(
					"Impossible to get info from QBlast service at this time. Check your network connection.\n");
		}
	}

	private URLConnection setQBlastServiceProperties(URLConnection conn) {

		URLConnection tmp = conn;

		conn.setDoOutput(true);
		conn.setUseCaches(false);

		tmp.setRequestProperty("User-Agent", "Biojava/RemoteQBlastService");
		tmp.setRequestProperty("Connection", "Keep-Alive");
		tmp.setRequestProperty("Content-type",
				"application/x-www-form-urlencoded");
		tmp.setRequestProperty("Content-length", "200");

		return tmp;
	}
    /** 
     * Set the tool identifier for QBlast. Defaults to 'biojavax'.
     * @param tool the new identifier.
     */
    public void setTool(String tool) {
        this.tool = tool;
    }

    /** 
     * Get the tool identifier for QBlast. Defaults to 'biojavax'.
     * @return the identifier.
     */
    public String getTool() {
        return this.tool;
    }

    /** 
     * Set the email for QBlast. Defaults to 'anonymous@biojava.org'.
     * @param email the new email.
     */
    public void setEmail(String email) {
        this.email = email;
    }

    /** 
     * Get the email for QBlast. Defaults to 'anonymous@biojava.org'.
     * @return the email.
     */
    public String getEmail() {
        return this.email;
    }	
}