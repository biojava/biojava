/*
 * BioJava development code
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
 * Author: Daniel Asarnow
 * Date:   2012-6-23
 */

package org.biojava.bio.structure.cath;

import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.io.util.FileDownloadUtils;
import org.biojava3.core.util.InputStreamProvider;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * @author Daniel Asarnow
 */
public class CathInstallation implements CathDatabase{

    public static final String DEFAULT_VERSION = "3.5.0";

    String cathVersion;

    public static final String domainListFileName = "CathDomainList";
    public static final String domainDescriptionFileName = "CathDomainDescriptionFile";
    public static final String nodeListFileName = "CathNames";
    public static final String domallFileName = "CathDomall";

    public static final String CATH_DOWNLOAD = "http://release.cathdb.info/";

    String cathDownloadUrl;

	public static final String NEWLINE;
	public static final String FILESPLIT ;

	static {
		NEWLINE     = System.getProperty("line.separator");
		FILESPLIT   = System.getProperty("file.separator");
	}

	String cacheLocation ;

    AtomicBoolean installedDomainList;
    AtomicBoolean installedDomainDescription;
    AtomicBoolean installedNodeList;
    AtomicBoolean installedDomall;

    final boolean useCathDomainDescriptionFile;
    final boolean parseCathFragments;

    Map<String, List<CathDomain>> pdbMap;
    Map<String, CathDomain> domainMap;
    Map<String, CathNode> cathTree;
    Map<String, List<CathFragment>> fragmentMap;
    
  
    	
    public CathInstallation(String cacheLocation, boolean usingCDDF, boolean parseCF) {
        setCacheLocation(cacheLocation);

        useCathDomainDescriptionFile = usingCDDF;
        parseCathFragments = parseCF;

        installedDomainDescription = new AtomicBoolean(false);
        installedDomainList = new AtomicBoolean(false);
        installedNodeList = new AtomicBoolean(false);
        installedDomall = new AtomicBoolean(false);

        cathVersion = DEFAULT_VERSION;
        cathDownloadUrl = CATH_DOWNLOAD;

        pdbMap = new HashMap<String, List<CathDomain>>();
        domainMap = new HashMap<String ,CathDomain>();
        cathTree = new HashMap<String, CathNode>();

        if (parseCathFragments) fragmentMap = new HashMap<String,List<CathFragment>>();

    }

    public CathInstallation(String cacheLocation) {
        this(cacheLocation, false, false);
    }

	public CathInstallation() {
		this((new UserConfiguration()).getPdbFilePath());
	}

    public String getDomainListFileName() {
        return cacheLocation + domainListFileName + ".v" + cathVersion;
    }

    public String getDomainDescriptionFileName() {
        return cacheLocation + domainDescriptionFileName + ".v" + cathVersion;
    }

    public String getNodeListFileName() {
        return cacheLocation + nodeListFileName + ".v" + cathVersion;
    }

    public String getDomallFileName() {
        return cacheLocation + domallFileName + ".v" + cathVersion;
    }

    public String getCathDownloadUrl() {
        return cathDownloadUrl;
    }

    public void setCathDownloadUrl(String cathDownloadUrl) {
        this.cathDownloadUrl = cathDownloadUrl;
    }

    public String getCacheLocation() {
        return cacheLocation;
    }

    public void setCacheLocation(String cacheLocation) {
        if ( !cacheLocation.endsWith(FILESPLIT) ) cacheLocation += FILESPLIT;
        this.cacheLocation = cacheLocation;
    }

    public AtomicBoolean getInstalledDomainList() {
        return installedDomainList;
    }

    public void setInstalledDomainList(AtomicBoolean installedDomainList) {
        this.installedDomainList = installedDomainList;
    }

    public AtomicBoolean getInstalledDomainDescription() {
        return installedDomainDescription;
    }

    public void setInstalledDomainDescription(AtomicBoolean installedDomainDescription) {
        this.installedDomainDescription = installedDomainDescription;
    }

    public AtomicBoolean getInstalledNodeList() {
        return installedNodeList;
    }

    public AtomicBoolean getInstalledDomall() {
        return installedDomall;
    }

    public void setInstalledNodeList(AtomicBoolean installedNodeList) {
        this.installedNodeList = installedNodeList;
    }

    public void setInstalledDomall(AtomicBoolean installedDomall) {
        this.installedDomall = installedDomall;
    }

    @Override
    public String getCathVersion() {
        return cathVersion;
    }

    @Override
    public CathNode getCathNode(String nodeId) {
        ensureNodeListInstalled();
        return cathTree.get(nodeId);
    }

    @Override
    public List<CathDomain> getByCategory(CathCategory category) {
        if (useCathDomainDescriptionFile) {
            ensureDomainDescriptionInstalled();
        } else {
            ensureDomallInstalled();
        }
        ensureNodeListInstalled();
        List<CathDomain> matches = new ArrayList<CathDomain>();
        CathNode node;
        for ( String nodeId : cathTree.keySet() ) {
            if ( (node = cathTree.get(nodeId)).getCategory() == category ) {
                matches.add( domainMap.get( node.getRepresentative() ) );
            }
        }
        return matches;
    }

    @Override
    public List<CathDomain> filterByCathCode(String query) {
        if (useCathDomainDescriptionFile) {
            ensureDomainDescriptionInstalled();
        } else {
            ensureDomallInstalled();
        }
        List<CathDomain> matches = new ArrayList<CathDomain>();
        for ( String k : domainMap.keySet() ) {
            if ( domainMap.get(k).getCATH().startsWith(query) ) {
                matches.add( domainMap.get(k) );
            }
        }
        return matches;
    }

    @Override
    public List<CathNode> getTree(CathDomain domain) {
        CathNode node = getCathNode( domain.getCATH() );
        List<CathNode> tree = new ArrayList<CathNode>();
        while (node != null) {
            node = getCathNode( node.getParentId() );
            if (node != null) tree.add(node);
        }
        Collections.reverse(tree);
        return tree;
    }

    @Override
    public List<CathDomain> filterByNodeName(String query) {
        ensureNodeListInstalled();
        List<CathNode> matchingNodes = new ArrayList<CathNode>();
        CathNode node;
        for ( String nodeId : cathTree.keySet() ) {
            if ( (node = cathTree.get(nodeId) ).getDescription().startsWith(query) ) {
                matchingNodes.add(node);
            }
        }
        List<CathDomain> matches = new ArrayList<CathDomain>();
        for (CathNode n : matchingNodes) {
            matches.addAll(getDomainsByNodeId(n.getNodeId()));
        }
        return matches;
    }

    @Override
    public List<CathDomain> filterByDescription(String query) {
        if (useCathDomainDescriptionFile) {
            ensureDomainDescriptionInstalled();
        } else {
            ensureDomallInstalled();
        }
        List<CathDomain> matches = new ArrayList<CathDomain>();
        for ( String k : domainMap.keySet() ) {
            if ( domainMap.get(k).getName().startsWith(query) ) {
                matches.add( domainMap.get(k) );
            }
        }
        return matches;
    }

    @Override
    public CathDomain getDescriptionByNodeId(String nodeId) {
        if (useCathDomainDescriptionFile) {
            ensureDomainDescriptionInstalled();
        } else {
            ensureDomallInstalled();
        }
        CathNode node = getCathNode(nodeId);
        return domainMap.get(node.getRepresentative());
    }

    @Override
    public List<CathDomain> getDomainsForPdb(String pdbId) {
        if (useCathDomainDescriptionFile) {
            ensureDomainDescriptionInstalled();
        } else {
            ensureDomallInstalled();
        }
        
      // cath IDs in lower case...
        return pdbMap.get(pdbId.toLowerCase());
    }

    @Override
    public CathDomain getDomainByCathId(String cathId) {
        if (useCathDomainDescriptionFile) {
            ensureDomainDescriptionInstalled();
        } else {
            ensureDomallInstalled();
        }
        return domainMap.get(cathId);
    }

    @Override
    public CathDomain getDescriptionByCathId(String cathId) {
        if (useCathDomainDescriptionFile) {
            ensureDomainDescriptionInstalled();
        } else {
            ensureDomallInstalled();
        }
        return domainMap.get(cathId);
    }

    @Override
    public List<CathDomain> getDomainsByNodeId(String nodeId) {
        if (useCathDomainDescriptionFile) {
            ensureDomainDescriptionInstalled();
        } else {
            ensureDomallInstalled();
        }
        List<CathDomain> domains = new ArrayList<CathDomain>();
        for (String domainName : domainMap.keySet()) {
            CathDomain description = domainMap.get(domainName);
            if ( description.getCATH().startsWith(nodeId) ) {
                domains.add(description);
            }
        }
        return domains;
    }

    @Override
    public List<CathFragment> getFragmentsByPdbId(String pdbId) {
        if ( useCathDomainDescriptionFile || !parseCathFragments ) return null;
        ensureDomallInstalled();
        return fragmentMap.get(pdbId);
    }

    private void parseCathDomainList() throws IOException {
		File file = new File(getDomainListFileName());
		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));
		parseCathDomainList(buffer);
    }

    private void parseCathDomainList(BufferedReader bufferedReader) throws IOException{
        String line;
     //   int counter = 0;
        while ( (line = bufferedReader.readLine()) != null ) {
            if ( line.startsWith("#") ) continue;
            CathDomain cathDomain = parseCathListFileLine(line);
           // counter++;
                        
            String pdbId = cathDomain.getPdbId().substring(0,4); // includes chain letter
			            
            List<CathDomain> domainList;
			if ( pdbMap.containsKey(pdbId)){
				domainList = pdbMap.get(pdbId);
			} else {
				domainList = new ArrayList<CathDomain>();
				pdbMap.put(pdbId,domainList);
			}

            domainList.add(cathDomain);

            domainMap.put( cathDomain.getDomainName(), cathDomain );
        }
    }

    private void parseCathNames() throws IOException {
		File file = new File(getNodeListFileName());
		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));
		parseCathNames(buffer);
    }

    private void parseCathNames(BufferedReader bufferedReader) throws IOException{
        String line;
        int counter = 0;
        while ( (line = bufferedReader.readLine()) != null ) {
            if ( line.startsWith("#") ) continue;
            CathNode cathNode = parseCathNamesFileLine(line);
            cathTree.put(cathNode.getNodeId(), cathNode);
        }
    }

    private void parseCathDomainDescriptionFile() throws IOException {
		File file = new File(getDomainDescriptionFileName());
		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));
		parseCathDomainDescriptionFile(buffer);
    }

    private void parseCathDomainDescriptionFile(BufferedReader bufferedReader) throws IOException{
        String line;
        DateFormat dateFormat = new SimpleDateFormat("dd-MMM-yyyy");
        int counter = 0;
        CathDomain cathDescription = null; //TODO initialize these or catch NPE
        StringBuilder name = null;
        StringBuilder source = null;
        StringBuilder seqh = null;
        StringBuilder seqs = null;
        List<CathSegment> segments = null;
        CathSegment segment = null;
        StringBuilder sseqh = null;
        StringBuilder sseqs = null;
        while ( (line = bufferedReader.readLine()) != null ) {        	
            if ( line.startsWith("#") ) continue;
            if ( line.startsWith("FORMAT") ) {
                cathDescription = new CathDomain();
                cathDescription.setFormat( line.substring(10) );

                name = new StringBuilder();
                source = new StringBuilder();
                seqh = new StringBuilder();
                seqs = new StringBuilder();

            } else if ( line.startsWith("DOMAIN") ) {
                cathDescription.setDomainName( line.substring(10) );
            } else if ( line.startsWith("VERSION") ) {
                cathDescription.setVersion( line.substring(10) );
            } else if ( line.startsWith("VERDATE") ) {
                try {
                    cathDescription.setDate( dateFormat.parse( line.substring(10) ) );
                } catch (ParseException e) {
                    e.printStackTrace();
                }
            } else if ( line.startsWith("NAME") ) {
                name.append( line.substring(10) );
            } else if ( line.startsWith("SOURCE") ) {
                source.append( line.substring(10) );
            } else if ( line.startsWith("CATHCODE") ) {
                cathDescription.setCATH( line.substring(10) );
            } else if ( line.startsWith("DLENGTH") ) {
                cathDescription.setLength( Integer.parseInt( line.substring(10) ) );
            } else if ( line.startsWith("DSEQH") ) {
                seqh.append( line.substring(10) );
            } else if ( line.startsWith("DSEQS") ) {
                seqs = seqs.append( line.substring(10) );
            } else if ( line.startsWith("NSEGMENTS") ) {
                segments = new ArrayList<CathSegment>();
            } else if ( line.startsWith("SEGMENT") ) {
                segment = new CathSegment();
                sseqh = new StringBuilder();
                sseqs = new StringBuilder();
            } else if ( line.startsWith("SRANGE") ) {
                int startStart = line.indexOf("=",10) + 1;
                int startStop = line.indexOf(" ",10);
                int stopStart = line.indexOf("=",startStop) + 1;
//                Integer start = Integer.parseInt( line.substring(startStart,startStop) );
//                Integer stop = Integer.parseInt( line.substring(stopStart, line.length()) );
                segment.setStart( line.substring(startStart,startStop) );
                segment.setStop( line.substring(stopStart) );
            } else if ( line.startsWith("SLENGTH") ) {
                segment.setLength( Integer.parseInt( line.substring(10) ) );
            } else if ( line.startsWith("SSEQH") ) {
                sseqh.append( line.substring(10) );
            } else if ( line.startsWith("SSEQS") ) {
                sseqs.append( line.substring(10) );
            } else if ( line.startsWith("ENDSEG") ) {
                segments.add( segment );
                segment.setSegmentId( segments.size() );
                segment.setSequenceHeader( sseqh.toString() );
                segment.setSequence( sseqs.toString() );
            } else if ( line.startsWith("//") ) {
                cathDescription.setName( name.toString() );
                cathDescription.setSource( source.toString() );
                cathDescription.setSequenceHeader( seqh.toString() );
                cathDescription.setSequence( seqs.toString() );
                cathDescription.setSegments(segments);
                counter++;

                String pdbId = cathDescription.getPdbId().substring(0,4); // includes chain letter
			    List<CathDomain> domainList;
			    if ( pdbMap.containsKey(pdbId)){
				    domainList = pdbMap.get(pdbId);
			    } else {
				    domainList = new ArrayList<CathDomain>();
				    pdbMap.put(pdbId,domainList);
			    }

                domainList.add(cathDescription);

                domainMap.put( cathDescription.getDomainName(), cathDescription );

            }
        }
//        transposeDomainData();
    }

/*    private void transposeDomainData() {
        ensureDomainListInstalled();
        for (String k : domainMap.keySet() ) {
            cathMap.get(k).getDomain().setResolution(domainMap.get(k).getResolution());
            cathMap.get(k).getDomain().setSOLID(domainMap.get(k).getSOILD());
        }
    }*/

    private CathDomain parseCathListFileLine(String line) {
        CathDomain cathDomain = new CathDomain();
        String [] token = line.split("\\s+");
        cathDomain.setDomainName(token[0]);
        cathDomain.setClassId(Integer.parseInt(token[1]));
        cathDomain.setArchitectureId(Integer.parseInt(token[2]));
        cathDomain.setTopologyId(Integer.parseInt(token[3]));
        cathDomain.setHomologyId(Integer.parseInt(token[4]));
        cathDomain.setSequenceFamilyId(Integer.parseInt(token[5]));
        cathDomain.setOrthologousSequenceFamilyId(Integer.parseInt(token[6]));
        cathDomain.setLikeSequenceFamilyId(Integer.parseInt(token[7]));
        cathDomain.setIdenticalSequenceFamilyId(Integer.parseInt(token[8]));
        cathDomain.setDomainCounter(Integer.parseInt(token[9]));
        cathDomain.setLength(Integer.parseInt(token[10]));
        cathDomain.setResolution(Double.parseDouble(token[11]));
        return cathDomain;
    }

    private CathNode parseCathNamesFileLine(String line) {
        CathNode cathNode = new CathNode();
        String[] token = line.split("\\s+",3);
        cathNode.setNodeId( token[0] );
        int idx = token[0].lastIndexOf(".");
        if ( idx == -1 ) idx = token[0].length();
        cathNode.setParentId( token[0].substring( 0, idx ) );
        cathNode.setRepresentative( token[1] );
        cathNode.setDescription( token[2].replace(":","") );
        return cathNode;
    }

    private void parseCathDomall() throws IOException{
		File file = new File(getDomallFileName());
		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));
		parseCathDomall(buffer);
    }

    private void parseCathDomall(BufferedReader bufferedReader) throws IOException{
        String line;
        while ( ((line = bufferedReader.readLine()) != null) ) {
            if ( line.startsWith("#") ) continue;
            if ( line.length() == 0 ) continue;
            String[] token = line.split("\\s+");
            String chainId = token[0];
            Integer numberOfDomains = Integer.parseInt( token[1].substring(1) );
            Integer numberOfFragments = Integer.parseInt( token[2].substring(1) );
            int domIdx = 3;
            int segIdx;
            Integer sstop;
            Integer sstart;
            Integer fstart;
            Integer fstop;
            Integer flength;
            for (int i=1; i<=numberOfDomains; i++) {
                DecimalFormat df = new DecimalFormat("00");
                String domainId;
                CathDomain domain;

//                This logic is necessary because singular domains may be labeled with 00 or 01.
//                If there is more than one domain, they are always numbered from 01.
                if (numberOfDomains==1) {
                    domainId = chainId + "00";
                    domain = domainMap.get(domainId);
                    if (domain==null) {
                        domainId = chainId + "01";
                        domain = domainMap.get(domainId);
                    }
                } else {
                    domainId = chainId + df.format(i);
                    domain = domainMap.get(domainId);
                }

                Integer numberOfSegments = Integer.parseInt( token[domIdx] );

                if ( domain == null ) {
                    domIdx += 6*numberOfSegments + 1;
                    continue;
                }

                List<CathSegment> segments = new ArrayList<CathSegment>(numberOfSegments);
                segIdx = 1; // Offset from domIdx.
                for (int j=1; j<=numberOfSegments; j++) {
                    CathSegment segment = new CathSegment();
                    segment.setSegmentId(j);
//                    String chainLetter = token[domIdx+segIdx]; // Redundant unless some domains cross chain boundaries.
                    sstart = Integer.parseInt( token[domIdx + segIdx + 1] );
                    String sstartInsertion = token[domIdx + segIdx + 2];
                    sstartInsertion = sstartInsertion.equals("-") ? "" : sstartInsertion;
//                    String chainLetter = token[domIdx+segIdx+4]; // Redundant unless some segments cross chain boundaries.
                    segment.setStart(sstart + sstartInsertion);

                    sstop = Integer.parseInt( token[domIdx + segIdx + 4] );
                    String sstopInsertion = token[domIdx + segIdx + 5];
                    sstopInsertion = sstopInsertion.equals("-") ? "" : sstopInsertion;

                    segment.setStart(sstart + sstartInsertion);
                    segment.setStop(sstop + sstopInsertion);
                    segment.setLength(1 + sstop - sstart);
                    segments.add(segment);

                    segIdx += 6;
                }
                domain.setSegments(segments);
                domIdx += 6*numberOfSegments + 1;
            }
            if (parseCathFragments) {
            List<CathFragment> fragments = new ArrayList<CathFragment>(numberOfFragments);
                for (int i=1; i<=numberOfFragments; i++) {
                    CathFragment fragment = new CathFragment();
                    fragment.setFragmentId(i);
//                    String chainLetter = token[domIdx]; // Redundant unless some fragments cross chain boundaries.
                    fstart = Integer.parseInt( token[domIdx+1] );
                    String fstartInsertion = token[domIdx + 2];
                    fstartInsertion = fstartInsertion.equals("-") ? "" : fstartInsertion;
                    fragment.setStart(fstart + fstartInsertion);
//                    String chainLetter = token[domIdx+3]; // Redundant unless some fragments cross chain boundaries.
                    fstop = Integer.parseInt( token[domIdx+4] );
                    String fstopInsertion = token[domIdx + 5];
                    fstopInsertion = fstopInsertion.equals("-") ? "" : fstopInsertion;
                    fragment.setStop(fstop + fstopInsertion);
                    flength = Integer.parseInt( token[domIdx + 6].replaceAll("[^0-9]","") );
                    fragment.setLength(flength);
                    fragments.add(fragment);
                    domIdx += 7;
                }
                fragmentMap.put(chainId, fragments);
            }
//            if ( domIdx != token.length ); // Problems.
        }
    }

    protected void downloadFileFromRemote(URL remoteURL, File localFile) throws FileNotFoundException, IOException{
        System.out.println("downloading " + remoteURL + " to: " + localFile);
        
        long timeS = System.currentTimeMillis();
    	File tempFile  = File.createTempFile(FileDownloadUtils.getFilePrefix(localFile), "."+ FileDownloadUtils.getFileExtension(localFile));
        
        FileOutputStream out = new FileOutputStream(tempFile);

        InputStream in = remoteURL.openStream();
        byte[] buf = new byte[4 * 1024]; // 4K buffer
        int bytesRead;
        while ((bytesRead = in.read(buf)) != -1) {
            out.write(buf, 0, bytesRead);
        }
        in.close();
        out.close();
        
        FileDownloadUtils.copy(tempFile,localFile);
        
        // delete the tmp file			
     	tempFile.delete();
        
     	long size =  localFile.length();
     	
     	double disp = size / 1024.0;
     	String unit = " kB";
     	if ( disp > 1024 ) {
     		unit = " MB";
     		disp = disp / 1024.0;
     	}
     	long timeE = System.currentTimeMillis();
        System.out.println("downloaded " + String.format("%.1f",disp) + unit  + " in " + (timeE - timeS)/1000 + " sec.");
	}

    private boolean domainDescriptionFileAvailable(){
		String fileName = getDomainDescriptionFileName();
		File f = new File(fileName);
		return f.exists();
	}

    private boolean domainListFileAvailable(){
		String fileName = getDomainListFileName();
		File f = new File(fileName);
		return f.exists();
	}

    private boolean nodeListFileAvailable(){
		String fileName = getNodeListFileName();
		File f = new File(fileName);
		return f.exists();
	}

    private boolean domallFileAvailable() {
        String fileName = getDomallFileName();
        File f= new File(fileName);
        return f.exists();
    }

	protected void downloadDomainListFile() throws FileNotFoundException, IOException{
		String remoteFilename = domainListFileName;
		URL url = new URL(cathDownloadUrl + "v" + cathVersion + "/" + remoteFilename);
		String localFileName = getDomainListFileName();
		File localFile = new File(localFileName);
		downloadFileFromRemote(url, localFile);
    }

	protected void downloadDomainDescriptionFile() throws FileNotFoundException, IOException{
		String remoteFilename = domainDescriptionFileName;
		URL url = new URL(cathDownloadUrl + "v" + cathVersion + "/" + remoteFilename);
		String localFileName = getDomainDescriptionFileName();
		File localFile = new File(localFileName);
		downloadFileFromRemote(url, localFile);
    }

	protected void downloadNodeListFile() throws FileNotFoundException, IOException{
		String remoteFilename = nodeListFileName;
		URL url = new URL(cathDownloadUrl + "v" + cathVersion + "/" + remoteFilename);
		String localFileName = getNodeListFileName();
		File localFile = new File(localFileName);
		downloadFileFromRemote(url, localFile);
    }

    protected void downloadDomallFile() throws IOException {
        String remoteFileName = domallFileName;
        URL url = new URL(cathDownloadUrl + "v" + cathVersion + "/" + remoteFileName);
        String localFileName = getDomallFileName();
        File localFile = new File(localFileName);
        downloadFileFromRemote(url, localFile);
    }

    public void ensureDomainListInstalled(){
		if ( installedDomainList.get() ) return;

		if ( ! domainListFileAvailable() ){
			try {
				downloadDomainListFile();
			} catch (Exception e){
				e.printStackTrace();
				installedDomainList.set(false);
				return;
			}
		}

		try {
			parseCathDomainList();
		} catch (Exception e){
			e.printStackTrace();
			installedDomainList.set(false);
			return;
		}
		installedDomainList.set(true);
	}

    public void ensureDomainDescriptionInstalled(){
		if ( installedDomainDescription.get() ) return;

		if ( ! domainDescriptionFileAvailable() ){
			try {
				downloadDomainDescriptionFile();
			} catch (Exception e){
				e.printStackTrace();
				installedDomainDescription.set(false);
				return;
			}
		}

		try {
			parseCathDomainDescriptionFile();
		} catch (Exception e){
			e.printStackTrace();
			installedDomainDescription.set(false);
			return;
		}
		installedDomainDescription.set(true);
	}

    public void ensureNodeListInstalled(){
		if ( installedNodeList.get() ) return;

		if ( ! nodeListFileAvailable() ){
			try {
				downloadNodeListFile();
			} catch (Exception e){
				e.printStackTrace();
				installedNodeList.set(false);
				return;
			}
		}

		try {
			parseCathNames();
		} catch (Exception e){
			e.printStackTrace();
			installedNodeList.set(false);
			return;
		}
		installedNodeList.set(true);
	}

    public void ensureDomallInstalled() {
        ensureDomainListInstalled();

        if ( !installedDomainList.get() ) {
            installedDomall.set(false);
            return;
        }

        if ( installedDomall.get() ) return;

        if ( ! domallFileAvailable() ){
            try {
                downloadDomallFile();
            } catch (Exception e) {
                e.printStackTrace();
                installedDomall.set(false);
                return;
            }
        }

        try {
            parseCathDomall();
        } catch (Exception e) {
            e.printStackTrace();
            installedDomall.set(false);
            return;
        }
        installedDomall.set(true);
    }

}
