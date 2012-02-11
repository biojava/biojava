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

import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.CMD;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.DATABASE;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.EMAIL;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.PROGRAM;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.QUERY;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.TOOL;
import static org.biojava3.ws.alignment.qblast.BlastOutputParameterEnum.RID;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.io.util.IOUtils;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.ws.alignment.RemotePairwiseAlignmentOutputProperties;
import org.biojava3.ws.alignment.RemotePairwiseAlignmentProperties;
import org.biojava3.ws.alignment.RemotePairwiseAlignmentService;

/**
 * Provides a simple way of submitting BLAST request to the QBlast service at NCBI.
 * <p>
 * NCBI provides a Blast server through a CGI-BIN interface. This service simply encapsulates an access to it by giving
 * users access to get/set methods to fix sequence, program and database as well as advanced options.
 * </p>
 * <p>
 * The philosophy behind this service is to disconnect submission of Blast requests from collection of Blast results.
 * This is done so to allow a user to submit multiple Blast requests while allowing recovery of the reports at a later
 * time.
 * </p>
 * <p>
 * Presently, only blastall programs are accessible.
 * </p>
 * 
 * @author Sylvain Foisy, Diploide BioIT
 * @author Gediminas Rimsa
 */
public class NCBIQBlastService implements RemotePairwiseAlignmentService {
	/**
	 * Number of milliseconds by which expected job execution time is incremented if it is not finished yet. Subsequent
	 * calls to {@link #isReady(String, long)} method will return false until at least this much time passes.
	 */
	public static final long WAIT_INCREMENT = 3000;

	private static final MapToStringTransformer MAP_TO_STRING_TRANSFORMER = new MapToStringTransformer();
	private static final String SERVICE_URL = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi";
	private static final String DEFAULT_EMAIL = "anonymous@biojava.org";
	private static final String DEFAULT_TOOL = "biojava3";

	private URL serviceUrl;
	private String email = DEFAULT_EMAIL;
	private String tool = DEFAULT_TOOL;

	private Map<String, BlastJob> jobs = new HashMap<String, BlastJob>();

	public NCBIQBlastService() {
		try {
			serviceUrl = new URL(SERVICE_URL);
		} catch (MalformedURLException e) {
			throw new RuntimeException("It looks like the URL for NCBI QBlast service is wrong. Cause: " + e.getMessage(), e);
		}
	}

	/**
	 * A simple method to check the availability of the QBlast service. Sends {@code Info} command to QBlast
	 * 
	 * @return QBlast info output concatenated to String
	 * @throws Exception if unable to connect to the NCBI QBlast service
	 */
	public String getRemoteBlastInfo() throws Exception {
		OutputStreamWriter writer = null;
		BufferedReader reader = null;
		try {
			URLConnection serviceConnection = setQBlastServiceProperties(serviceUrl.openConnection());
			writer = new OutputStreamWriter(serviceConnection.getOutputStream());
			writer.write("CMD=Info");
			writer.flush();
			reader = new BufferedReader(new InputStreamReader(serviceConnection.getInputStream()));
			StringBuilder sb = new StringBuilder();
			String line;
			while ((line = reader.readLine()) != null) {
				sb.append(line);
				sb.append(System.getProperty("line.separator"));
			}
			return sb.toString();
		} catch (IOException e) {
			throw new Exception("Impossible to get info from QBlast service at this time. Cause: " + e.getMessage(), e);
		} finally {
			IOUtils.close(reader);
			IOUtils.close(writer);
		}
	}

	/**
	 * Converts given sequence to String and calls
	 * {@link #sendAlignmentRequest(String, RemotePairwiseAlignmentProperties)}
	 */
	@Override
	public String sendAlignmentRequest(Sequence<Compound> seq, RemotePairwiseAlignmentProperties rpa) throws Exception {
		return sendAlignmentRequest(seq.getSequenceAsString(), rpa);
	}

	/**
	 * Converts given GenBank GID to String and calls
	 * {@link #sendAlignmentRequest(String, RemotePairwiseAlignmentProperties)}
	 */
	public String sendAlignmentRequest(int gid, RemotePairwiseAlignmentProperties rpa) throws Exception {
		return sendAlignmentRequest(Integer.toString(gid), rpa);
	}

	/**
	 * Sends the Blast request via the Put command of the CGI-BIN interface. Uses all of the parameters specified in
	 * {@code alignmentProperties} (parameters PROGRAM and DATABASE are required).
	 * 
	 * @param query : a {@code String} representing a sequence or Genbank ID
	 * @param alignmentProperties : a {@code RemotePairwiseAlignmentProperties} object representing alignment properties
	 * @return the request id for this sequence, necessary to fetch results after completion
	 * @throws Exception if unable to connect to the NCBI QBlast service or if no sequence or required parameters
	 *             PROGRAM and DATABASE are not set
	 */
	@Override
	public String sendAlignmentRequest(String query, RemotePairwiseAlignmentProperties alignmentProperties) throws Exception {
		Map<String, String> params = new HashMap<String, String>();
		for (String key : alignmentProperties.getAlignmentOptions()) {
			params.put(key, alignmentProperties.getAlignmentOption(key));
		}

		if (query == null || query.isEmpty()) {
			throw new IllegalArgumentException("Impossible to execute QBlast request. The sequence has not been set.");
		}
		if (!params.containsKey(PROGRAM.name())) {
			throw new IllegalArgumentException("Impossible to execute QBlast request. Parameter PROGRAM has not been set.");
		}
		if (!params.containsKey(DATABASE.name())) {
			throw new IllegalArgumentException("Impossible to execute QBlast request. Parameter DATABASE has not been set.");
		}

		params.put(CMD.name(), "Put");
		params.put(QUERY.name(), query);
		params.put(TOOL.name(), getTool());
		params.put(EMAIL.name(), getEmail());

		String putCmd = MAP_TO_STRING_TRANSFORMER.transform(params);

		OutputStreamWriter writer = null;
		BufferedReader reader = null;
		try {
			URLConnection serviceConnection = setQBlastServiceProperties(serviceUrl.openConnection());
			writer = new OutputStreamWriter(serviceConnection.getOutputStream());
			writer.write(putCmd);
			writer.flush();

			// Get the response
			reader = new BufferedReader(new InputStreamReader(serviceConnection.getInputStream()));

			// find request id and time of execution
			BlastJob job = new BlastJob();
			String line;
			while ((line = reader.readLine()) != null) {
				if (!line.contains("class=\"error\"") && !line.contains("Message ID#")) {
					// if there is no error, capture RID and RTOE
					if (line.contains("RID = ")) {
						String[] arr = line.split("=");
						job.setId(arr[1].trim());
					} else if (line.contains("RTOE = ")) {
						String[] arr = line.split("=");
						job.setStartTimestamp(System.currentTimeMillis());
						job.setExpectedExecutionTime(Long.parseLong(arr[1].trim()) * 1000);
					}
					jobs.put(job.getId(), job);
				} else {
					// handle QBlast error message

					// Capture everything to the left of this HTML statement...
					String[] tmp = line.split("</p></li></ul>");

					// Only the error message is on the right side of this...
					String[] moreTmp = tmp[0].split("<p class=\"error\">");
					throw new Exception("NCBI QBlast refused this request because: " + moreTmp[1].trim());
				}

			}
			if (job != null && job.getId() != null) {
				return job.getId();
			} else {
				throw new Exception("Unable to retrieve request ID");
			}
		} catch (IOException e) {
			throw new IOException("An error occured submiting sequence to BLAST server. Cause: " + e.getMessage(), e);
		} finally {
			IOUtils.close(reader);
			IOUtils.close(writer);
		}
	}

	/**
	 * Wrapper method for {@link #isReady(String, long)}, omitting unnecessary {@code present} property.
	 * 
	 * @see #isReady(String, long)
	 */
	public boolean isReady(String id) throws Exception {
		return isReady(id, 0);
	}

	/**
	 * Checks for completion of request.
	 * <p/>
	 * If expected execution time (RTOE) is available for request, this method will always return false until that time
	 * passes. This is done to prevent sending unnecessary requests to the server.
	 * 
	 * @param id : request id, which was returned by {@code sendAlignmentRequest} method
	 * @param present : is not used, can be any value
	 * @return a boolean value telling if the request has been completed
	 */
	@Override
	public boolean isReady(String id, long present) throws Exception {
		BlastJob job = jobs.get(id);
		if (job != null) {
			long expectedJobFinishTime = job.getStartTimestamp() + job.getExpectedExecutionTime();
			if (System.currentTimeMillis() < expectedJobFinishTime) {
				return false;
			}
		} else {
			// it might be a valid job from another session
			job = new BlastJob();
			job.setId(id);
			job.setStartTimestamp(System.currentTimeMillis());
			job.setExpectedExecutionTime(0);
		}

		OutputStreamWriter writer = null;
		BufferedReader reader = null;
		try {
			String checkRequest = "CMD=Get&RID=" + job.getId();
			URLConnection serviceConnection = setQBlastServiceProperties(serviceUrl.openConnection());
			writer = new OutputStreamWriter(serviceConnection.getOutputStream());
			writer.write(checkRequest);
			writer.flush();
			reader = new BufferedReader(new InputStreamReader(serviceConnection.getInputStream()));

			String line;
			while ((line = reader.readLine()) != null) {
				if (line.contains("READY")) {
					jobs.put(job.getId(), job);
					return true;
				} else if (line.contains("WAITING")) {
					job.setExpectedExecutionTime(job.getExpectedExecutionTime() + WAIT_INCREMENT);
					jobs.put(job.getId(), job);
					return false;
				} else if (line.contains("UNKNOWN")) {
					throw new IllegalArgumentException("Unknown request id - no results exist for it. Given id = " + id);
				}
			}
			return false;
		} catch (IOException ioe) {
			throw new Exception("It is not possible to fetch Blast report from NCBI at this time. Cause: " + ioe.getMessage(), ioe);
		} finally {
			IOUtils.close(reader);
			IOUtils.close(writer);
		}
	}

	/**
	 * Extracts the actual Blast report for given request id according to options provided in {@code outputProperties}
	 * argument.
	 * <p/>
	 * If the results are not ready yet, sleeps until they are available. If sleeping is not desired, call this method
	 * after {@code isReady} returns true
	 * 
	 * @param id : request id, which was returned by {@code sendAlignmentRequest} method
	 * @param outputProperties : an object specifying output formatting options
	 * @return an {@code InputStream} of results
	 * @throws Exception if it is not possible to recover the results
	 */
	public InputStream getAlignmentResults(String id, RemotePairwiseAlignmentOutputProperties outputProperties) throws Exception {
		Map<String, String> params = new HashMap<String, String>();
		for (String key : outputProperties.getOutputOptions()) {
			params.put(key, outputProperties.getOutputOption(key));
		}
		OutputStreamWriter writer = null;

		while (!isReady(id)) {
			Thread.sleep(WAIT_INCREMENT + 100);
		}

		params.put(CMD.name(), "Get");
		params.put(RID.name(), id);
		params.put(TOOL.name(), getTool());
		params.put(EMAIL.name(), getEmail());
		String getCmd = MAP_TO_STRING_TRANSFORMER.transform(params);

		try {
			URLConnection serviceConnection = setQBlastServiceProperties(serviceUrl.openConnection());
			writer = new OutputStreamWriter(serviceConnection.getOutputStream());
			writer.write(getCmd);
			writer.flush();
			return serviceConnection.getInputStream();
		} catch (IOException ioe) {
			throw new Exception("It is not possible to fetch Blast report from NCBI at this time. Cause: " + ioe.getMessage(), ioe);
		} finally {
			IOUtils.close(writer);
		}
	}

	/**
	 * Sends a delete request for given request id. Optional operation, ignores IOExceptions.<br/>
	 * Can be used after results of given search are no longer needed to be kept on Blast server
	 * 
	 * @param id request id, as returned by {@code sendAlignmentRequest} method
	 */
	public void sendDeleteRequest(String id) {
		OutputStreamWriter writer = null;
		try {
			String deleteRequest = "CMD=Delete&RID=" + id;
			URLConnection serviceConnection = setQBlastServiceProperties(serviceUrl.openConnection());
			writer = new OutputStreamWriter(serviceConnection.getOutputStream());
			writer.write(deleteRequest);
			writer.flush();
		} catch (IOException ignore) {
			// ignore it this is an optional operation
		} finally {
			IOUtils.close(writer);
		}
	}

	/**
	 * Sets properties for given URLConnection
	 * 
	 * @param conn URLConnection to set properties for
	 * @return given object after setting properties
	 */
	private URLConnection setQBlastServiceProperties(URLConnection conn) {
		conn.setDoOutput(true);
		conn.setUseCaches(false);
		conn.setRequestProperty("User-Agent", "Biojava/NCBIQBlastService");
		conn.setRequestProperty("Connection", "Keep-Alive");
		conn.setRequestProperty("Content-type", "application/x-www-form-urlencoded");
		conn.setRequestProperty("Content-length", "2000");
		return conn;
	}

	/**
	 * Set the tool identifier for QBlast. Defaults to {@value #DEFAULT_TOOL}
	 * 
	 * @param tool the new identifier
	 */
	public void setTool(String tool) {
		this.tool = tool;
	}

	/**
	 * Get the tool identifier for QBlast. Defaults to {@value #DEFAULT_TOOL}
	 * 
	 * @return the identifier
	 */
	public String getTool() {
		return this.tool;
	}

	/**
	 * Set the email for QBlast. Defaults to {@value #DEFAULT_EMAIL}
	 * 
	 * @param email the new email
	 */
	public void setEmail(String email) {
		this.email = email;
	}

	/**
	 * Get the email for QBlast. Defaults to {@value #DEFAULT_EMAIL}.
	 * 
	 * @return the email
	 */
	public String getEmail() {
		return this.email;
	}
}