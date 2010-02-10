// $Id: GoTerm.java,v 1.11 2009/11/10 23:12:22 cmzmasek Exp $
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

package org.forester.go;

import java.util.List;

import org.forester.phylogeny.data.PhylogenyData;

public interface GoTerm extends PhylogenyData, Comparable<GoTerm> {

    public List<GoId> getAltIds();

    public String getComment();

    public String getDefinition();

    public GoId getGoId();

    public GoNameSpace getGoNameSpace();

    public List<GoRelationship> getGoRelationships();

    public List<GoSubset> getGoSubsets();

    public List<GoXRef> getGoXRefs();

    public String getName();

    public List<GoId> getSuperGoIds();

    public boolean isObsolete();
}
