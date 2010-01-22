// $Id: HmmPfamOutputParser.java,v 1.16 2009/11/10 19:57:08 cmzmasek Exp $
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.io.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.surfacing.BasicDomain;
import org.forester.surfacing.BasicProtein;
import org.forester.surfacing.Domain;
import org.forester.surfacing.DomainId;
import org.forester.surfacing.Protein;
import org.forester.surfacing.SurfacingUtil;
import org.forester.util.ForesterUtil;

public final class HmmPfamOutputParser {

    private static final String     RETRO                       = "RETRO";
    private static final String     PHAGE                       = "PHAGE";
    private static final String     VIR                         = "VIR";
    private static final String     TRANSPOS                    = "TRANSPOS";
    private static final String     RV                          = "RV";
    private static final String     GAG                         = "GAG_";
    private static final String     HCV                         = "HCV_";                                                    // New. Added on Jun 11, after 1st submission.
    private static final String     HERPES                      = "Herpes_";                                                 // New. Added on Jun 11, after 1st submission.
    private static final int        E_VALUE_MAXIMUM_DEFAULT     = -1;
    private static final ReturnType RETURN_TYPE_DEFAULT         = ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN;
    private static final boolean    IGNORE_DUFS_DEFAULT         = false;
    private static final int        MAX_ALLOWED_OVERLAP_DEFAULT = -1;
    private final Set<DomainId>     _filter;
    private final FilterType        _filter_type;
    private final File              _input_file;
    private final String            _species;
    private final String            _model_type;
    private double                  _e_value_maximum;
    private Map<String, String>     _individual_domain_score_cutoffs;
    private boolean                 _ignore_dufs;
    private boolean                 _ignore_virus_like_ids;
    private boolean                 _allow_non_unique_query;
    private boolean                 _verbose;
    private int                     _max_allowed_overlap;
    private boolean                 _ignore_engulfed_domains;
    private ReturnType              _return_type;
    private int                     _proteins_encountered;
    private int                     _proteins_ignored_due_to_filter;
    private int                     _proteins_stored;
    private int                     _domains_encountered;
    private int                     _domains_ignored_due_to_duf;
    private int                     _domains_ignored_due_to_overlap;
    private int                     _domains_ignored_due_to_e_value;
    private int                     _domains_ignored_due_to_individual_score_cutoff;
    private int                     _domains_stored;
    private SortedSet<DomainId>     _domains_stored_set;
    private long                    _time;
    private int                     _domains_ignored_due_to_negative_domain_filter;
    private Map<String, Integer>    _domains_ignored_due_to_negative_domain_filter_counts_map;
    private int                     _domains_ignored_due_to_virus_like_id;
    private Map<String, Integer>    _domains_ignored_due_to_virus_like_id_counts_map;

    public HmmPfamOutputParser( final File input_file, final String species, final String model_type ) {
        _input_file = input_file;
        _species = species;
        _model_type = model_type;
        _filter = null;
        _filter_type = FilterType.NONE;
        init();
    }

    public HmmPfamOutputParser( final File input_file,
                                final String species,
                                final String model_type,
                                final Set<DomainId> filter,
                                final FilterType filter_type ) {
        _input_file = input_file;
        _species = species;
        _model_type = model_type;
        _filter = filter;
        _filter_type = filter_type;
        init();
    }

    private void actuallyAddProtein( final List<Protein> proteins, final Protein current_protein ) {
        final List<Domain> l = current_protein.getProteinDomains();
        for( final Domain d : l ) {
            getDomainsStoredSet().add( d.getDomainId() );
        }
        proteins.add( current_protein );
        ++_proteins_stored;
    }

    private void addProtein( final List<Protein> proteins, final Protein current_protein ) {
        if ( ( getFilterType() == FilterType.POSITIVE_PROTEIN ) || ( getFilterType() == FilterType.NEGATIVE_PROTEIN ) ) {
            final Set<DomainId> domain_ids_in_protein = new HashSet<DomainId>();
            for( final Domain d : current_protein.getProteinDomains() ) {
                domain_ids_in_protein.add( d.getDomainId() );
            }
            domain_ids_in_protein.retainAll( getFilter() );
            if ( getFilterType() == FilterType.POSITIVE_PROTEIN ) {
                if ( domain_ids_in_protein.size() > 0 ) {
                    actuallyAddProtein( proteins, current_protein );
                }
                else {
                    ++_proteins_ignored_due_to_filter;
                }
            }
            else {
                if ( domain_ids_in_protein.size() < 1 ) {
                    actuallyAddProtein( proteins, current_protein );
                }
                else {
                    ++_proteins_ignored_due_to_filter;
                }
            }
        }
        else {
            actuallyAddProtein( proteins, current_protein );
        }
    }

    public int getDomainsEncountered() {
        return _domains_encountered;
    }

    public int getDomainsIgnoredDueToDuf() {
        return _domains_ignored_due_to_duf;
    }

    public int getDomainsIgnoredDueToEval() {
        return _domains_ignored_due_to_e_value;
    }

    public int getDomainsIgnoredDueToIndividualScoreCutoff() {
        return _domains_ignored_due_to_individual_score_cutoff;
    }

    public int getDomainsIgnoredDueToNegativeDomainFilter() {
        return _domains_ignored_due_to_negative_domain_filter;
    }

    public Map<String, Integer> getDomainsIgnoredDueToNegativeDomainFilterCountsMap() {
        return _domains_ignored_due_to_negative_domain_filter_counts_map;
    }

    public int getDomainsIgnoredDueToOverlap() {
        return _domains_ignored_due_to_overlap;
    }

    public Map<String, Integer> getDomainsIgnoredDueToVirusLikeIdCountsMap() {
        return _domains_ignored_due_to_virus_like_id_counts_map;
    }

    public int getDomainsIgnoredDueToVirusLikeIds() {
        return _domains_ignored_due_to_virus_like_id;
    }

    public int getDomainsStored() {
        return _domains_stored;
    }

    public SortedSet<DomainId> getDomainsStoredSet() {
        return _domains_stored_set;
    }

    private double getEValueMaximum() {
        return _e_value_maximum;
    }

    private Set<DomainId> getFilter() {
        return _filter;
    }

    private FilterType getFilterType() {
        return _filter_type;
    }

    private Map<String, String> getIndividualDomainScoreCutoffs() {
        return _individual_domain_score_cutoffs;
    }

    private File getInputFile() {
        return _input_file;
    }

    private int getMaxAllowedOverlap() {
        return _max_allowed_overlap;
    }

    private String getModelType() {
        return _model_type;
    }

    public int getProteinsEncountered() {
        return _proteins_encountered;
    }

    public int getProteinsIgnoredDueToFilter() {
        return _proteins_ignored_due_to_filter;
    }

    public int getProteinsStored() {
        return _proteins_stored;
    }

    private ReturnType getReturnType() {
        return _return_type;
    }

    private String getSpecies() {
        return _species;
    }

    public long getTime() {
        return _time;
    }

    private void init() {
        _e_value_maximum = HmmPfamOutputParser.E_VALUE_MAXIMUM_DEFAULT;
        setIgnoreDufs( HmmPfamOutputParser.IGNORE_DUFS_DEFAULT );
        setReturnType( HmmPfamOutputParser.RETURN_TYPE_DEFAULT );
        _max_allowed_overlap = HmmPfamOutputParser.MAX_ALLOWED_OVERLAP_DEFAULT;
        setIndividualDomainScoreCutoffs( null );
        setIgnoreEngulfedDomains( false );
        setIgnoreVirusLikeIds( false );
        setAllowNonUniqueQuery( false );
        setVerbose( false );
        intitCounts();
    }

    private void intitCounts() {
        setDomainsStoredSet( new TreeSet<DomainId>() );
        setDomainsEncountered( 0 );
        setProteinsEncountered( 0 );
        setProteinsIgnoredDueToFilter( 0 );
        setDomainsIgnoredDueToNegativeFilter( 0 );
        setDomainsIgnoredDueToDuf( 0 );
        setDomainsIgnoredDueToEval( 0 );
        setDomainsIgnoredDueToIndividualScoreCutoff( 0 );
        setDomainsIgnoredDueToVirusLikeId( 0 );
        setDomainsIgnoredDueToOverlap( 0 );
        setDomainsStored( 0 );
        setProteinsStored( 0 );
        setTime( 0 );
        setDomainsIgnoredDueToVirusLikeIdCountsMap( new TreeMap<String, Integer>() );
        setDomainsIgnoredDueToNegativeDomainFilterCountsMap( new TreeMap<String, Integer>() );
    }

    private boolean isAllowNonUniqueQuery() {
        return _allow_non_unique_query;
    }

    private boolean isIgnoreDufs() {
        return _ignore_dufs;
    }

    private boolean isIgnoreEngulfedDomains() {
        return _ignore_engulfed_domains;
    }

    private boolean isIgnoreVirusLikeIds() {
        return _ignore_virus_like_ids;
    }

    private boolean isVerbose() {
        return _verbose;
    }

    public List<Protein> parse() throws IOException {
        intitCounts();
        final Set<String> queries = new HashSet<String>();
        final String error = ForesterUtil.isReadableFile( getInputFile() );
        if ( !ForesterUtil.isEmpty( error ) ) {
            throw new IOException( error );
        }
        final BufferedReader br = new BufferedReader( new FileReader( getInputFile() ) );
        String line;
        final List<Protein> proteins = new ArrayList<Protein>();
        Protein current_protein = null;
        int line_number = 0;
        boolean saw_double_slash = true;
        boolean can_parse_domains = false;
        boolean saw_parsed_for_domains = false;
        boolean saw_query_sequence = false;
        boolean was_not_unique = false;
        final long start_time = new Date().getTime();
        while ( ( line = br.readLine() ) != null ) {
            line_number++;
            if ( line.length() < 1 ) {
                continue;
            }
            else if ( line.startsWith( "Query sequence:" ) ) {
                ++_proteins_encountered;
                if ( !saw_double_slash ) {
                    throw new IOException( "unexpected format [line " + line_number + "] in ["
                            + getInputFile().getCanonicalPath() + "]" );
                }
                saw_double_slash = false;
                saw_query_sequence = true;
                was_not_unique = false;
                final String query = line.substring( 16 ).trim();
                if ( ForesterUtil.isEmpty( query ) ) {
                    throw new IOException( "query sequence cannot be empty [line " + line_number + "] in ["
                            + getInputFile().getCanonicalPath() + "]" );
                }
                if ( queries.contains( query ) ) {
                    if ( !isAllowNonUniqueQuery() ) {
                        throw new IOException( "query \"" + query + "\" is not unique [line " + line_number + "] in ["
                                + getInputFile().getCanonicalPath() + "]" );
                    }
                    else if ( isVerbose() ) {
                        ForesterUtil.printWarningMessage( getClass().getName(), "query \"" + query
                                + "\" is not unique [line " + line_number + "] in ["
                                + getInputFile().getCanonicalPath() + "]" );
                    }
                }
                else {
                    queries.add( query );
                }
                if ( current_protein != null ) {
                    throw new IOException( "unexpected format [line " + line_number + "] in ["
                            + getInputFile().getCanonicalPath() + "]" );
                }
                if ( getReturnType() == ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN ) {
                    current_protein = new BasicProtein( query, getSpecies() );
                }
                else {
                    throw new IllegalArgumentException( "unknown return type" );
                }
            }
            else if ( line.startsWith( "Accession:" ) ) {
                if ( !saw_query_sequence || ( current_protein == null ) ) {
                    throw new IOException( "unexpected format [line " + line_number + "] in ["
                            + getInputFile().getCanonicalPath() + "]" );
                }
                ( ( BasicProtein ) current_protein ).setAccession( line.substring( 11 ).trim() );
            }
            else if ( line.startsWith( "Description:" ) ) {
                if ( !saw_query_sequence || ( current_protein == null ) ) {
                    throw new IOException( "unexpected format [line " + line_number + "] in ["
                            + getInputFile().getCanonicalPath() + "]" );
                }
                if ( was_not_unique ) {
                    if ( getReturnType() == ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN ) {
                        current_protein = new BasicProtein( current_protein.getProteinId() + " "
                                + line.substring( 13 ).trim(), getSpecies() );
                    }
                }
                else {
                    ( ( BasicProtein ) current_protein ).setDescription( line.substring( 13 ).trim() );
                }
            }
            else if ( line.startsWith( "Parsed for domains:" ) ) {
                if ( !saw_query_sequence ) {
                    throw new IOException( "unexpected format [line " + line_number + "] in ["
                            + getInputFile().getCanonicalPath() + "]" );
                }
                saw_query_sequence = false;
                saw_parsed_for_domains = true;
            }
            else if ( saw_parsed_for_domains && line.startsWith( "--------" ) ) {
                can_parse_domains = true;
                saw_parsed_for_domains = false;
            }
            else if ( line.startsWith( "Alignments of top-scoring domains:" ) ) {
                if ( !can_parse_domains ) {
                    throw new IOException( "unexpected format [line " + line_number + "] in ["
                            + getInputFile().getCanonicalPath() + "]" );
                }
                can_parse_domains = false;
            }
            else if ( line.startsWith( "//" ) ) {
                can_parse_domains = false;
                saw_double_slash = true;
                if ( current_protein.getProteinDomains().size() > 0 ) {
                    if ( ( getMaxAllowedOverlap() != HmmPfamOutputParser.MAX_ALLOWED_OVERLAP_DEFAULT )
                            || isIgnoreEngulfedDomains() ) {
                        final int domains_count = current_protein.getNumberOfProteinDomains();
                        current_protein = SurfacingUtil.removeOverlappingDomains( getMaxAllowedOverlap(),
                                                                                  isIgnoreEngulfedDomains(),
                                                                                  current_protein );
                        final int domains_removed = domains_count - current_protein.getNumberOfProteinDomains();
                        _domains_stored -= domains_removed;
                        _domains_ignored_due_to_overlap += domains_removed;
                    }
                    addProtein( proteins, current_protein );
                }
                current_protein = null;
            }
            else if ( can_parse_domains && ( line.indexOf( "[no hits above thresholds]" ) == -1 ) ) {
                final String[] s = line.split( "\\s+" );
                if ( s.length != 10 ) {
                    throw new IOException( "unexpected format in hmmpfam output:  \"" + line + "\" [line "
                            + line_number + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                final String id = s[ 0 ];
                final String domain_count_str = s[ 1 ];
                final String from_str = s[ 2 ];
                final String to_str = s[ 3 ];
                final String query_match_str = s[ 4 ];
                final String hmm_match_str = s[ 7 ];
                final String score_str = s[ 8 ];
                final String e_value_str = s[ 9 ];
                int from = -1;
                int to = -1;
                double e_value = -1;
                double score = -1;
                boolean is_complete_hmm_match = false;
                boolean is_complete_query_match = false;
                try {
                    from = Integer.valueOf( from_str ).intValue();
                }
                catch ( final NumberFormatException e ) {
                    throw new IOException( "could not parse seq-f from \"" + line + "\" [line " + line_number
                            + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                try {
                    to = Integer.valueOf( to_str ).intValue();
                }
                catch ( final NumberFormatException e ) {
                    throw new IOException( "could not parse seq-t from \"" + line + "\" [line " + line_number
                            + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                try {
                    score = Double.valueOf( score_str ).doubleValue();
                }
                catch ( final NumberFormatException e ) {
                    throw new IOException( "could not parse score from \"" + line + "\" [line " + line_number
                            + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                try {
                    e_value = Double.valueOf( e_value_str ).doubleValue();
                }
                catch ( final NumberFormatException e ) {
                    throw new IOException( "could not parse E-value from \"" + line + "\" [line " + line_number
                            + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                if ( hmm_match_str.equals( "[]" ) ) {
                    is_complete_hmm_match = true;
                }
                else if ( !( hmm_match_str.equals( ".]" ) || hmm_match_str.equals( "[." ) || hmm_match_str
                        .equals( ".." ) ) ) {
                    throw new IOException( "unexpected format in hmmpfam output:  \"" + line + "\" [line "
                            + line_number + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                if ( query_match_str.equals( ".." ) ) {
                    is_complete_query_match = true;
                }
                else if ( !( query_match_str.equals( ".]" ) || query_match_str.equals( "[." ) || query_match_str
                        .equals( "[]" ) ) ) {
                    throw new IOException( "unexpected format in hmmpfam output:  \"" + line + "\" [line "
                            + line_number + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                final String[] ct = domain_count_str.split( "/" );
                if ( ct.length != 2 ) {
                    throw new IOException( "unexpected format in hmmpfam output:  \"" + line + "\" [line "
                            + line_number + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                final String number_str = ct[ 0 ];
                final String total_str = ct[ 1 ];
                int number = -1;
                int total = -1;
                try {
                    number = Integer.valueOf( ( number_str ) ).intValue();
                }
                catch ( final NumberFormatException e ) {
                    throw new IOException( "could not parse domain number from \"" + line + "\" [line " + line_number
                            + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                try {
                    total = Integer.valueOf( ( total_str ) ).intValue();
                }
                catch ( final NumberFormatException e ) {
                    throw new IOException( "could not parse domain count from \"" + line + "\" [line " + line_number
                            + "] in [" + getInputFile().getCanonicalPath() + "]" );
                }
                ++_domains_encountered;
                boolean failed_cutoff = false;
                if ( getIndividualDomainScoreCutoffs() != null ) {
                    if ( getIndividualDomainScoreCutoffs().containsKey( id ) ) {
                        final double cutoff = Double.parseDouble( getIndividualDomainScoreCutoffs().get( id ) );
                        if ( score < cutoff ) {
                            failed_cutoff = true;
                        }
                    }
                    else {
                        throw new IOException( "could not find a score cutoff value for domain id \"" + id
                                + "\" [line " + line_number + "] in [" + getInputFile().getCanonicalPath() + "]" );
                    }
                }
                final String uc_id = id.toUpperCase();
                if ( failed_cutoff ) {
                    ++_domains_ignored_due_to_individual_score_cutoff;
                }
                else if ( ( getEValueMaximum() != HmmPfamOutputParser.E_VALUE_MAXIMUM_DEFAULT )
                        && ( e_value > getEValueMaximum() ) ) {
                    ++_domains_ignored_due_to_e_value;
                }
                else if ( isIgnoreDufs() && uc_id.startsWith( "DUF" ) ) {
                    ++_domains_ignored_due_to_duf;
                }
                else if ( isIgnoreVirusLikeIds()
                        && ( uc_id.contains( VIR ) || uc_id.contains( PHAGE ) || uc_id.contains( RETRO )
                                || uc_id.contains( TRANSPOS ) || uc_id.startsWith( RV ) || uc_id.startsWith( GAG )
                                || uc_id.startsWith( HCV ) || uc_id.startsWith( HERPES ) ) ) {
                    ForesterUtil.increaseCountingMap( getDomainsIgnoredDueToVirusLikeIdCountsMap(), id );
                    ++_domains_ignored_due_to_virus_like_id;
                }
                else if ( ( getFilterType() == FilterType.NEGATIVE_DOMAIN )
                        && getFilter().contains( new DomainId( id ) ) ) {
                    ++_domains_ignored_due_to_negative_domain_filter;
                    ForesterUtil.increaseCountingMap( getDomainsIgnoredDueToNegativeDomainFilterCountsMap(), id );
                }
                else {
                    final BasicDomain pd = new BasicDomain( id,
                                                            from,
                                                            to,
                                                            ( short ) number,
                                                            ( short ) total,
                                                            e_value,
                                                            score );
                    current_protein.addProteinDomain( pd );
                    ++_domains_stored;
                }
            }
        } // while ( ( line = br.readLine() ) != null )
        setTime( new Date().getTime() - start_time );
        if ( !saw_double_slash ) {
            throw new IOException( "file ends unexpectedly [line " + line_number + "]" );
        }
        return proteins;
    }

    public void setAllowNonUniqueQuery( final boolean allow_non_unique_query ) {
        _allow_non_unique_query = allow_non_unique_query;
    }

    private void setDomainsEncountered( final int domains_encountered ) {
        _domains_encountered = domains_encountered;
    }

    private void setDomainsIgnoredDueToDuf( final int domains_ignored_due_to_duf ) {
        _domains_ignored_due_to_duf = domains_ignored_due_to_duf;
    }

    public void setDomainsIgnoredDueToEval( final int domains_ignored_due_to_e_value ) {
        _domains_ignored_due_to_e_value = domains_ignored_due_to_e_value;
    }

    public void setDomainsIgnoredDueToIndividualScoreCutoff( final int domains_ignored_due_to_individual_score_cutoff ) {
        _domains_ignored_due_to_individual_score_cutoff = domains_ignored_due_to_individual_score_cutoff;
    }

    private void setDomainsIgnoredDueToNegativeDomainFilterCountsMap( final Map<String, Integer> domains_ignored_due_to_negative_domain_filter_counts_map ) {
        _domains_ignored_due_to_negative_domain_filter_counts_map = domains_ignored_due_to_negative_domain_filter_counts_map;
    }

    private void setDomainsIgnoredDueToNegativeFilter( final int domains_ignored_due_to_negative_domain_filter ) {
        _domains_ignored_due_to_negative_domain_filter = domains_ignored_due_to_negative_domain_filter;
    }

    private void setDomainsIgnoredDueToOverlap( final int domains_ignored_due_to_overlap ) {
        _domains_ignored_due_to_overlap = domains_ignored_due_to_overlap;
    }

    private void setDomainsIgnoredDueToVirusLikeId( final int i ) {
        _domains_ignored_due_to_virus_like_id = i;
    }

    private void setDomainsIgnoredDueToVirusLikeIdCountsMap( final Map<String, Integer> domains_ignored_due_to_virus_like_id_counts_map ) {
        _domains_ignored_due_to_virus_like_id_counts_map = domains_ignored_due_to_virus_like_id_counts_map;
    }

    private void setDomainsStored( final int domains_stored ) {
        _domains_stored = domains_stored;
    }

    private void setDomainsStoredSet( final SortedSet<DomainId> _storeddomains_stored ) {
        _domains_stored_set = _storeddomains_stored;
    }

    public void setEValueMaximum( final double e_value_maximum ) {
        if ( e_value_maximum < 0.0 ) {
            throw new IllegalArgumentException( "attempt to set the maximum E-value to a negative value" );
        }
        _e_value_maximum = e_value_maximum;
    }

    public void setIgnoreDufs( final boolean ignore_dufs ) {
        _ignore_dufs = ignore_dufs;
    }

    /**
     * To ignore domains which are completely engulfed by domains (individual
     * ones or stretches of overlapping ones) with better support values.
     * 
     * 
     * @param ignored_engulfed_domains
     */
    public void setIgnoreEngulfedDomains( final boolean ignore_engulfed_domains ) {
        _ignore_engulfed_domains = ignore_engulfed_domains;
    }

    public void setIgnoreVirusLikeIds( final boolean ignore_virus_like_ids ) {
        _ignore_virus_like_ids = ignore_virus_like_ids;
    }

    /**
     * Sets the individual domain score cutoff values (for example, gathering
     * thresholds from Pfam). Domain ids are the keys, cutoffs the values.
     * 
     * @param individual_domain_score_cutoffs
     */
    public void setIndividualDomainScoreCutoffs( final Map<String, String> individual_domain_score_cutoffs ) {
        _individual_domain_score_cutoffs = individual_domain_score_cutoffs;
    }

    public void setMaxAllowedOverlap( final int max_allowed_overlap ) {
        if ( max_allowed_overlap < 0 ) {
            throw new IllegalArgumentException( "Attempt to set max allowed overlap to less than zero." );
        }
        _max_allowed_overlap = max_allowed_overlap;
    }

    private void setProteinsEncountered( final int proteins_encountered ) {
        _proteins_encountered = proteins_encountered;
    }

    private void setProteinsIgnoredDueToFilter( final int proteins_ignored_due_to_filter ) {
        _proteins_ignored_due_to_filter = proteins_ignored_due_to_filter;
    }

    private void setProteinsStored( final int proteins_stored ) {
        _proteins_stored = proteins_stored;
    }

    public void setReturnType( final ReturnType return_type ) {
        _return_type = return_type;
    }

    private void setTime( final long time ) {
        _time = time;
    }

    public void setVerbose( final boolean verbose ) {
        _verbose = verbose;
    }

    public static enum FilterType {
        NONE, POSITIVE_PROTEIN, NEGATIVE_PROTEIN, NEGATIVE_DOMAIN
    }

    public static enum ReturnType {
        UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN
    }
}
