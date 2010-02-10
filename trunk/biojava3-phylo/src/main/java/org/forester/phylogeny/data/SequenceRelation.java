// $Id: SequenceRelation.java,v 1.16 2009/11/17 03:51:34 cmzmasek Exp $
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

public class SequenceRelation implements PhylogenyData {

    @Override
    public StringBuffer asSimpleText() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public StringBuffer asText() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public PhylogenyData copy() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        // TODO Auto-generated method stub
    }

    public static enum SEQUENCE_RELATION_TYPE {
        orthology, one_to_one_orthology, super_orthology, paralogy, ultra_paralogy, xenology, unknown, other;
    }
}
