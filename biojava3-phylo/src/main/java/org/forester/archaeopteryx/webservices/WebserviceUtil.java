// $Id: WebserviceUtil.java,v 1.4 2009/12/12 00:14:39 cmzmasek Exp $
// forester -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2010 Christian M. Zmasek
// Copyright (C) 2008-2010 Burnham Institute for Medical Research
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

package org.forester.archaeopteryx.webservices;

import java.util.ArrayList;
import java.util.List;

import org.forester.archaeopteryx.webservices.WebservicesManager.WsPhylogenyFormat;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterUtil;
import org.forester.util.ForesterUtil.PhylogenyNodeField;

public final class WebserviceUtil {

    public static final String TAX_CODE_TO_SCI_NAME = "tax_code_to_sci_name";
    public static final String TREE_FAM_INST        = "tree_fam";
    public static final String PFAM_INST            = "pfam";
    public static final String TOL_WEBSERVER        = "http://tolweb.org/onlinecontributors/app?service=external&page=xml/TreeStructureService&node_id="
                                                            + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER;
    public static final String TOL_NAME             = "Tree of Life";
    public static final String TREE_BASE_NAME       = "TreeBase";
    public static final String TREE_FAM_NAME        = "TreeFam";
    public static final String PFAM_NAME            = "Pfam";
    public static final String PFAM_SERVER          = "http://pfam.janelia.org";

    public static List<PhylogeniesWebserviceClient> createDefaultClients() {
        final List<PhylogeniesWebserviceClient> clients = new ArrayList<PhylogeniesWebserviceClient>();
        clients
                .add( new BasicPhylogeniesWebserviceClient( TOL_NAME,
                                                            "Read Tree from Tree of Life...",
                                                            "Use ToL webservice to obtain a phylogeny",
                                                            "Please enter a Tree of Life node identifier\n(Examples: "
                                                                    + "19386 for Cephalopoda, 2461 for Cnidaria, 2466 for Deuterostomia)",
                                                            WsPhylogenyFormat.TOL_XML_RESPONSE,
                                                            PhylogenyNodeField.TAXONOMY_SCIENTIFIC_NAME,
                                                            WebserviceUtil.TOL_WEBSERVER,
                                                            true,
                                                            "http://tolweb.org",
                                                            null ) );
        clients
                .add( new BasicPhylogeniesWebserviceClient( TREE_BASE_NAME,
                                                            "Read Tree from TreeBase...",
                                                            "Use TreeBase to obtain a phylogeny",
                                                            "Please enter a TreeBase tree identifier\n(Examples: 2654, 825, 3306, 2518, 2406)",
                                                            WsPhylogenyFormat.NH,
                                                            PhylogenyNodeField.TAXONOMY_SCIENTIFIC_NAME,
                                                            "http://130.132.27.194/treebase/TreeBASE.acgi?PickedItems=Tree"
                                                                    + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER
                                                                    + "&Button=ATV.nhx",
                                                            true,
                                                            "http://www.treebase.org",
                                                            TAX_CODE_TO_SCI_NAME ) );
        clients
                .add( new BasicPhylogeniesWebserviceClient( PFAM_NAME,
                                                            "Read Gene Tree from Pfam...",
                                                            "Use  Pfam to obtain a (full) gene tree",
                                                            "Please enter a Pfam (PF) accession number\n(Examples: 01849 for NAC, 00452 for Bcl-2, 00046 for Homeobox)",
                                                            WsPhylogenyFormat.PFAM,
                                                            null,
                                                            PFAM_SERVER + "/family/tree/download?alnType=full&acc=PF"
                                                                    + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER,
                                                            false,
                                                            PFAM_SERVER,
                                                            PFAM_INST ) );
        clients
                .add( new BasicPhylogeniesWebserviceClient( TREE_FAM_NAME,
                                                            "Read Full Gene Tree from TreeFam...",
                                                            "Use TreeFam to obtain a (full) gene tree",
                                                            "Please enter a TreeFam (TF) accession number\n(Examples: 101004 for Cyclin D, 315938 for Hox, 105310 for Wnt)",
                                                            WsPhylogenyFormat.NHX,
                                                            null,
                                                            "http://www.treefam.org/cgi-bin/getdata.pl?ac=TF"
                                                                    + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER
                                                                    + "&f=full.nhx",
                                                            true,
                                                            "http://www.treefam.org",
                                                            TREE_FAM_INST ) );
        clients
                .add( new BasicPhylogeniesWebserviceClient( TREE_FAM_NAME,
                                                            "Read Clean Gene Tree from TreeFam...",
                                                            "Use TreeFam to obtain a (\"clean\") gene tree",
                                                            "Please enter a TreeFam (TF) accession number\n(Examples: 101004 for Cyclin D, 315938 for Hox, 105310 for Wnt)",
                                                            WsPhylogenyFormat.NHX,
                                                            null,
                                                            "http://www.treefam.org/cgi-bin/getdata.pl?ac=TF"
                                                                    + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER
                                                                    + "&f=clean.nhx",
                                                            true,
                                                            "http://www.treefam.org",
                                                            TREE_FAM_INST ) );
        return clients;
    }

    static void extractSpTremblAccFromNodeName( final Phylogeny phy, final String source ) {
        final PreorderTreeIterator it = new PreorderTreeIterator( phy );
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( !ForesterUtil.isEmpty( n.getNodeName() ) ) {
                final String name = n.getNodeName();
                final int i = name.lastIndexOf( "/" );
                if ( i > 0 ) {
                    final String acc_str = name.substring( 0, i );
                    if ( !ForesterUtil.isEmpty( acc_str ) ) {
                        final Sequence seq = new Sequence();
                        final Accession acc = new Accession( acc_str, source );
                        seq.setAccession( acc );
                        n.getNodeData().setSequence( seq );
                    }
                }
            }
        }
    }

    public static void processInstructions( final PhylogeniesWebserviceClient client, final Phylogeny phylogeny ) {
        if ( client.getProcessingInstructions().equals( WebserviceUtil.TAX_CODE_TO_SCI_NAME ) ) {
            WebserviceUtil.transferTaxonomyCodeToScientificName( phylogeny );
        }
        else if ( client.getProcessingInstructions().equals( WebserviceUtil.TREE_FAM_INST ) ) {
            WebserviceUtil.transferInternalTaxonomyCodeToScientificName( phylogeny );
            WebserviceUtil.transferExternalScientificNameToTaxonomyCode( phylogeny );
            WebserviceUtil.transferSequenceNameToSequenceAccession( phylogeny, "ensembl" );
            WebserviceUtil.setTaxonomyIdentifierType( phylogeny, "ncbi" );
        }
        else if ( client.getProcessingInstructions().equals( WebserviceUtil.PFAM_INST ) ) {
            WebserviceUtil.extractSpTremblAccFromNodeName( phylogeny, "sptrembl" );
        }
    }

    static void setTaxonomyIdentifierType( final Phylogeny phy, final String type ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() && ( n.getNodeData().getTaxonomy().getIdentifier() != null ) ) {
                n.getNodeData().getTaxonomy().setIdentifier( new Identifier( n.getNodeData().getTaxonomy()
                        .getIdentifier().getValue(), type ) );
            }
        }
    }

    static void transferExternalScientificNameToTaxonomyCode( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() && n.getNodeData().isHasTaxonomy() ) {
                final String name = n.getNodeData().getTaxonomy().getScientificName();
                if ( !ForesterUtil.isEmpty( name ) && PhyloXmlUtil.TAXOMONY_CODE_PATTERN.matcher( name ).matches() ) {
                    n.getNodeData().getTaxonomy().setScientificName( "" );
                    n.getNodeData().getTaxonomy().setTaxonomyCode( name );
                }
            }
        }
    }

    static void transferInternalTaxonomyCodeToScientificName( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && n.getNodeData().isHasTaxonomy() ) {
                final String name = n.getNodeData().getTaxonomy().getTaxonomyCode();
                if ( !ForesterUtil.isEmpty( name ) ) {
                    n.getNodeData().getTaxonomy().setScientificName( name );
                    n.getNodeData().getTaxonomy().setTaxonomyCode( "" );
                }
            }
        }
    }

    static void transferSequenceNameToSequenceAccession( final Phylogeny phy, final String source ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasSequence() ) {
                final String name = n.getNodeData().getSequence().getName();
                if ( !ForesterUtil.isEmpty( name ) ) {
                    n.getNodeData().getSequence().setName( "" );
                    n.getNodeData().getSequence().setAccession( new Accession( name, source ) );
                }
            }
        }
    }

    static void transferTaxonomyCodeToScientificName( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                final String name = n.getNodeData().getTaxonomy().getTaxonomyCode();
                if ( !ForesterUtil.isEmpty( name ) ) {
                    n.getNodeData().getTaxonomy().setScientificName( name );
                    n.getNodeData().getTaxonomy().setTaxonomyCode( "" );
                }
            }
        }
    }
}
