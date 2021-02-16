/* This file is part of PyMesh. Copyright (c) 2016 by Qingnan Zhou */
#pragma once
#ifndef GEOGRAMBASE_H
#define GEOGRAMBASE_H

#include <memory>

#include <geogram/basic/common.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>

using namespace GEO;

using GeoMeshPtr = std::shared_ptr<Mesh>;

/**
 * Base class for all GeoGram wrappers.
 * It takes care of initializing GeoGram environment.
 */
class GeogramBase {
    public:
        GeogramBase() {
            initialize();
            CmdLine::import_arg_group("standard");
            CmdLine::import_arg_group("pre");
            CmdLine::import_arg_group("remesh");
            CmdLine::import_arg_group("algo");
            CmdLine::import_arg_group("post");
            CmdLine::import_arg_group("opt");
            CmdLine::import_arg_group("co3ne");
            CmdLine::import_arg_group("tet");
            CmdLine::import_arg_group("poly");
            Logger::instance()->set_pretty(true);
        }

        virtual ~GeogramBase() = default;
};

#endif