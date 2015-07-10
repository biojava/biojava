package org.biojava.nbio.ws.uniprot;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.http.HttpEntity;
import org.apache.http.HttpStatus;
import org.apache.http.NameValuePair;
import org.apache.http.client.entity.UrlEncodedFormEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.client.methods.HttpRequestBase;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.message.BasicNameValuePair;
import org.apache.http.util.EntityUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import static org.biojava.nbio.ws.uniprot.UniprotConstants.*;

class BatchRetrieveUniprotEntries implements BatchRetrieve {
	
	private final static Logger logger = LoggerFactory.getLogger(BatchRetrieveUniprotEntries.class);
	
	@Override
	public void retrieve(Iterable<String> ids, BatchRetrievalFormat format, ResponseConsumer out) throws IOException {
		JobIdConsumer consumer = new JobIdConsumer(out);
		String jobId = createJob(ids, consumer).getJobId();
		if (!INVALID_JOB_ID.equals(jobId))
			retrieve(jobId, format, out);
	}

	private JobIdConsumer createJob(Iterable<String> ids, JobIdConsumer out) throws IOException {
		makeNetworkRequest(prepareHttpPost(ids), out);
		return out;
	}
	
	private void makeNetworkRequest(HttpRequestBase base, ResponseConsumer out) throws IOException{
		CloseableHttpClient httpclient = HttpClients.createDefault();
		CloseableHttpResponse response = httpclient.execute(base);
		try {
			int status = response.getStatusLine().getStatusCode();
			HttpEntity entity = response.getEntity();
			out.setReader(new InputStreamReader(response.getEntity().getContent()));
			if (HttpStatus.SC_OK == status) {
				out.onSucess();
			} else {
				out.onFail(status);
			}
		    EntityUtils.consume(entity);
		} finally {
		    response.close();
		}
	}
	
	private void retrieve(String jobId, BatchRetrievalFormat format, ResponseConsumer out) throws IOException {
		String queryURL = String.format((BatchRetrievalFormat.fastaWithIsoforms == format ? queryUrlWithIsoforms : queryUrl), jobId, format.getExtenstion());
		makeNetworkRequest(new HttpGet(queryURL), out);
	}

	private HttpPost prepareHttpPost(Iterable<String> ids) throws UnsupportedEncodingException {
		HttpPost httpPost = new HttpPost(uploadListURL);
		List <NameValuePair> nvps = new ArrayList <NameValuePair>();
		nvps.add(new BasicNameValuePair("format", "job"));
		nvps.add(new BasicNameValuePair("from", "ACC"));
		nvps.add(new BasicNameValuePair("to", "ACC"));
		nvps.add(new BasicNameValuePair("landingPage", "false"));
		StringBuilder builder = prepareIdString(ids);
		nvps.add(new BasicNameValuePair("uploadQuery", builder.toString()));
		httpPost.setEntity(new UrlEncodedFormEntity(nvps));
		return httpPost;
	}

	private StringBuilder prepareIdString(Iterable<String> ids) {
		StringBuilder builder = new StringBuilder();
		Iterator<String> iterator = ids.iterator();
		while(iterator.hasNext()) {
			builder.append(String.format("%s%s", iterator.next(), (iterator.hasNext() ? " " : "")));
		}
		return builder;
	}

	
	
	/**
	 * Utility to cread the stream into a 
	 * string
	 * @param reader
	 * @return the string that is read or an empty string
	 */
	public static String readAsString(Reader reader) throws Exception
	{
		StringBuffer buffer = new StringBuffer();
		char[] buf = new char[1024];
		int numRead = 0;
		while ((numRead = reader.read(buf)) != -1)
		{
			buffer.append(buf, 0, numRead);
		}
		return buffer.toString();
	}
	
	
	/**
	 * Implementation of NetworkResponseConsumer to consume a job id 
	 * responded from the server
	 * @author pbansal
	 *
	 */
	private static class JobIdConsumer implements ResponseConsumer {
		private String jobId = null;
		private Reader reader = null;
		private ResponseConsumer consumer;
		
		public JobIdConsumer(ResponseConsumer consumer) {
			this.consumer = consumer;
		}
		
		@Override
		public void onSucess() {
			try {
				jobId = readAsString(reader);
			} catch(Exception e) {
				logger.error("Error while retrieving the jobId ", e.getMessage());
			}
		}

		@Override
		public void onFail(int httpStatusCode) {
			consumer.onFail(httpStatusCode);
		}

		@Override
		public void setReader(Reader reader) {
			this.reader = reader;
		}
		
		public String getJobId() {
			return jobId == null ? INVALID_JOB_ID : jobId;
		}
	}
	
} // BatchRetrieveUniprotEntries
