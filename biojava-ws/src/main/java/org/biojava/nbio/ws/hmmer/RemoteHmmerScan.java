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
package org.biojava.nbio.ws.hmmer;

import net.sf.json.JSONArray;
import net.sf.json.JSONObject;
import org.biojava.nbio.core.sequence.ProteinSequence;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLEncoder;
import java.util.SortedSet;
import java.util.TreeSet;


/** Makes remote calls to the Hmmer3 web site and returns Pfam domain annotations for an input protein sequence.
 * 
 * @author Andreas Prlic
 * @since 3.0.3
 */
public class RemoteHmmerScan implements HmmerScan {

	public static String HMMER_SERVICE = "http://hmmer.janelia.org/search/hmmscan";

	public boolean DEFAULT_SEARCH_CUT_GA = true;

	private boolean searchWithCutGA;
	public RemoteHmmerScan(){
		searchWithCutGA = DEFAULT_SEARCH_CUT_GA;
	}

	@Override
	public  SortedSet<HmmerResult> scan(ProteinSequence sequence) throws IOException {

		URL url = new URL(HMMER_SERVICE);

		return scan(sequence, url);

	}

	/** Scans a protein sequence for Pfam profile matches.
	 * 
	 * @param sequence
	 * @param serviceLocation
	 * @return
	 * @throws IOException
	 */
	public SortedSet<HmmerResult> scan(ProteinSequence sequence, URL serviceLocation) throws IOException{


		String urlparameters = prepareParameters(sequence, searchWithCutGA);

		HttpURLConnection connection = (HttpURLConnection) serviceLocation.openConnection();
		connection.setDoOutput(true);
		connection.setDoInput(true);
		connection.setConnectTimeout(15000); // 15 sec
		connection.setInstanceFollowRedirects(false);
		connection.setRequestMethod("POST");
		connection.setRequestProperty("Content-Type", "application/x-www-form-urlencoded");
		connection.setRequestProperty("Accept", "application/json");

		connection.setRequestProperty("Content-Length", "" +
				Integer.toString(urlparameters.toString().getBytes().length));

		//Send request
		DataOutputStream wr = new DataOutputStream (
				connection.getOutputStream ());
		wr.writeBytes (urlparameters);
		wr.flush ();
		wr.close ();

		//Now get the redirect URL
		URL respUrl = new URL( connection.getHeaderField( "Location" ));

		int responseCode = connection.getResponseCode();
		if ( responseCode == 500){
			System.err.println("something went wrong!" + urlparameters);
			System.err.println(connection.getResponseMessage());
		}

		HttpURLConnection connection2 = (HttpURLConnection) respUrl.openConnection();
		connection2.setRequestMethod("GET");
		connection2.setRequestProperty("Accept", "application/json");


		//Get the response 
		BufferedReader in = new BufferedReader(
				new InputStreamReader(
						connection2.getInputStream()));

		String inputLine;

		StringBuffer result = new StringBuffer();
		while ((inputLine = in.readLine()) != null)
			result.append(inputLine);
		in.close();

		// process the response and build up a container for the data.

		SortedSet<HmmerResult> results = new TreeSet<HmmerResult>();
		try {
			JSONObject json =  JSONObject.fromObject(result.toString());

			JSONObject hmresults = json.getJSONObject("results");


			JSONArray hits = hmresults.getJSONArray("hits");

			for(int i =0 ; i < hits.size() ; i++){
				JSONObject hit = hits.getJSONObject(i);
				//System.out.println("hit: "+ hit);

				HmmerResult hmmResult = new HmmerResult();

				Object dclO = hit.get("dcl");
				Integer dcl = -1;
				if ( dclO instanceof Long){
					Long dclL = (Long) dclO;
					dcl = dclL.intValue();
				} else if ( dclO instanceof Integer){
					dcl = (Integer) dclO;
				} 
				
				hmmResult.setAcc((String)hit.get("acc"));
				hmmResult.setDcl(dcl);
				hmmResult.setDesc((String)hit.get("desc"));
				hmmResult.setEvalue(Float.parseFloat((String)hit.get("evalue")));
				hmmResult.setName((String)hit.get("name"));
				hmmResult.setNdom((Integer)hit.get("ndom"));
				hmmResult.setNreported((Integer)hit.get("nreported"));
				hmmResult.setPvalue((Double)hit.get("pvalue"));
				hmmResult.setScore(Float.parseFloat((String)hit.get("score")));

				JSONArray hmmdomains = hit.getJSONArray("domains");

				SortedSet<HmmerDomain> domains = new TreeSet<HmmerDomain>();
				for ( int j= 0 ; j < hmmdomains.size() ; j++){
					JSONObject d = hmmdomains.getJSONObject(j);
					//System.out.println(d);
					Integer is_reported = getInteger(d.get("is_reported"));
					if ( is_reported != 1) {
						//System.out.println("excluding: " + d);
						continue;
					}

					HmmerDomain dom = new HmmerDomain();
					dom.setAliLenth((Integer)d.get("aliL"));
					dom.setHmmAcc((String)d.get("alihmmacc"));
					dom.setHmmDesc((String)d.get("alihmmdesc"));

					dom.setHmmFrom(getInteger(d.get("alihmmfrom")));
					dom.setHmmTo(getInteger(d.get("alihmmto")));
					dom.setSimCount((Integer)d.get("aliSimCount"));
					dom.setSqFrom(getInteger(d.get("alisqfrom")));
					dom.setSqTo(getInteger(d.get("alisqto")));
					dom.setHmmName((String)d.get("alihmmname"));
					domains.add(dom);
				}

				hmmResult.setDomains(domains);

				results.add(hmmResult);
			}
		} catch (Exception e){
			e.printStackTrace();
		}

		return results;

	}

	private String prepareParameters(ProteinSequence sequence,
			boolean searchWithCutGA) throws UnsupportedEncodingException  {
		StringBuffer urlParameters = new StringBuffer();



		urlParameters.append("hmmdb=");
		urlParameters.append(URLEncoder.encode("pfam","UTF-8"));


		// search by cut_ga, not e-values
		//cut_ga' : '',
		if ( searchWithCutGA ) {
			urlParameters.append("&cut_ga=");
			//urlParameters.append(URLEncoder.encode("1","UTF-8"));


		}

		urlParameters.append("&seq=");

		urlParameters.append(URLEncoder.encode(">seq\n"+sequence.toString(), "UTF-8"));




		return urlParameters.toString();
	}

	private Integer getInteger(Object object) {
		if ( object instanceof Integer)
			return (Integer) object;
		else if( object instanceof String)
			return Integer.parseInt((String) object);

		return null;
	}

}
