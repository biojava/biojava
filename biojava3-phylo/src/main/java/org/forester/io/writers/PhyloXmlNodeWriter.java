// $Id: PhyloXmlNodeWriter.java,v 1.21 2009/12/12 00:14:39 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2000-2009 Christian M. Zmasek
// Copyright (C) 2007-2009 Burnham Institute for Medical Research
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

package org.forester.io.writers;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.util.ForesterUtil;

public class PhyloXmlNodeWriter {

    public static void toPhyloXml( final Writer w, final PhylogenyNode node, final int level, final String indentation )
            throws IOException {
        String ind = "";
        if ( indentation.length() > 0 ) {
            ind = indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE;
        }
        if ( !ForesterUtil.isEmpty( node.getNodeName() ) ) {
            PhylogenyDataUtil.appendElement( w, PhyloXmlMapping.NODE_NAME, node.getNodeName(), indentation );
        }
        if ( node.getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) {
            PhylogenyDataUtil.appendElement( w, PhyloXmlMapping.BRANCH_LENGTH, String.valueOf( ForesterUtil.round( node
                    .getDistanceToParent(), PhyloXmlUtil.ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT ) ), indentation );
        }
        if ( node.getBranchData() != null ) {
            node.getBranchData().toPhyloXML( w, level, ind );
        }
        if ( node.getNodeData() != null ) {
            node.getNodeData().toPhyloXML( w, level, ind );
        }
    }
}
