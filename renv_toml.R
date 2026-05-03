#!/usr/bin/env Rscript
# TOML-based dependency management for R projects
# Provides pyproject.toml-like functionality for R

# Install required packages if not available
if (!requireNamespace("RcppTOML", quietly = TRUE)) {
  install.packages("RcppTOML", repos = "https://cloud.r-project.org")
}

library(RcppTOML)

#' Read and parse Rproject.toml
#' @param path Path to Rproject.toml file
#' @return Parsed TOML configuration
read_rproject_toml <- function(path = "Rproject.toml") {
  if (!file.exists(path)) {
    stop("Rproject.toml not found. Run: create_rproject_toml()")
  }
  
  config <- parseTOML(path)
  return(config)
}

#' Install dependencies from Rproject.toml
#' @param groups Character vector of dependency groups to install
#' @param config Optional pre-parsed TOML config
install_from_toml <- function(groups = c("dependencies"), config = NULL) {
  if (is.null(config)) {
    config <- read_rproject_toml()
  }
  
  cat("Installing dependencies from Rproject.toml...\n")
  
  # Ensure renv is available
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv", repos = "https://cloud.r-project.org")
  }
  
  all_packages <- list()
  
  for (group in groups) {
    if (group %in% names(config)) {
      deps <- config[[group]]
      cat(sprintf("Installing %s group:\n", group))
      
      for (pkg_name in names(deps)) {
        pkg_spec <- deps[[pkg_name]]
        
        # Handle different package specification formats
        if (is.list(pkg_spec)) {
          # Complex specification: { source = "bioconductor", version = "1.0" }
          source <- pkg_spec$source %||% "cran"
          version <- pkg_spec$version %||% ""
          
          if (source == "bioconductor") {
            package_ref <- paste0("bioc::", pkg_name)
          } else if (source == "github") {
            package_ref <- pkg_spec$repo %||% pkg_name
          } else {
            package_ref <- pkg_name
          }
          
          if (version != "" && version != "*") {
            package_ref <- paste0(package_ref, "@", version)
          }
        } else {
          # Simple specification: "package" = "*" or "package" = "1.0.0"
          if (pkg_spec == "*" || pkg_spec == "") {
            package_ref <- pkg_name
          } else {
            package_ref <- paste0(pkg_name, "@", pkg_spec)
          }
        }
        
        cat(sprintf("  Installing: %s\n", package_ref))
        
        tryCatch({
          renv::install(package_ref)
          all_packages[[pkg_name]] <- package_ref
        }, error = function(e) {
          cat(sprintf("  Warning: Failed to install %s: %s\n", package_ref, e$message))
        })
      }
    }
  }
  
  # Install GitHub dependencies if present
  if ("github-dependencies" %in% names(config)) {
    cat("Installing GitHub dependencies:\n")
    github_deps <- config[["github-dependencies"]]
    
    for (repo in names(github_deps)) {
      ref <- github_deps[[repo]]
      package_ref <- if (ref == "main" || ref == "master") repo else paste0(repo, "@", ref)
      
      cat(sprintf("  Installing: %s\n", package_ref))
      tryCatch({
        renv::install(package_ref)
        all_packages[[basename(repo)]] <- package_ref
      }, error = function(e) {
        cat(sprintf("  Warning: Failed to install %s: %s\n", package_ref, e$message))
      })
    }
  }
  
  cat("\nCreating renv snapshot...\n")
  renv::snapshot(prompt = FALSE)
  
  cat(sprintf("Successfully installed %d packages.\n", length(all_packages)))
  return(invisible(all_packages))
}

#' Add a new dependency to Rproject.toml
#' @param package Package name
#' @param version Version constraint (default: "*")
#' @param group Dependency group (default: "dependencies")
#' @param source Package source (default: "cran")
add_dependency <- function(package, version = "*", group = "dependencies", source = "cran") {
  config <- read_rproject_toml()
  
  # Ensure the group exists
  if (!group %in% names(config)) {
    config[[group]] <- list()
  }
  
  # Add the package
  if (source == "cran" && version == "*") {
    config[[group]][[package]] <- "*"
  } else {
    config[[group]][[package]] <- list(
      version = version,
      source = source
    )
  }
  
  # Write back to file
  writeLines(toTOML(config), "Rproject.toml")
  
  cat(sprintf("Added %s to %s group in Rproject.toml\n", package, group))
  
  # Install the package
  cat("Installing package...\n")
  package_ref <- if (source == "bioconductor") {
    paste0("bioc::", package)
  } else {
    package
  }
  
  if (version != "*") {
    package_ref <- paste0(package_ref, "@", version)
  }
  
  renv::install(package_ref)
  renv::snapshot(prompt = FALSE)
  
  cat("Package installed and lockfile updated.\n")
}

#' Remove a dependency from Rproject.toml
#' @param package Package name
#' @param group Dependency group (default: "dependencies")
remove_dependency <- function(package, group = "dependencies") {
  config <- read_rproject_toml()
  
  if (group %in% names(config) && package %in% names(config[[group]])) {
    config[[group]][[package]] <- NULL
    
    # Write back to file
    writeLines(toTOML(config), "Rproject.toml")
    
    cat(sprintf("Removed %s from %s group in Rproject.toml\n", package, group))
    
    # Remove from renv
    renv::remove(package)
    renv::snapshot(prompt = FALSE)
    
    cat("Package removed and lockfile updated.\n")
  } else {
    cat(sprintf("Package %s not found in %s group.\n", package, group))
  }
}

#' List dependencies from Rproject.toml
#' @param group Specific group to list (default: all groups)
list_dependencies <- function(group = NULL) {
  config <- read_rproject_toml()
  
  dep_groups <- c("dependencies", "dev-dependencies", "optional-dependencies", "github-dependencies")
  
  if (!is.null(group)) {
    dep_groups <- intersect(dep_groups, group)
  }
  
  for (grp in dep_groups) {
    if (grp %in% names(config)) {
      cat(sprintf("\n[%s]\n", grp))
      deps <- config[[grp]]
      
      for (pkg in names(deps)) {
        spec <- deps[[pkg]]
        if (is.list(spec)) {
          cat(sprintf("  %s: %s (source: %s)\n", 
                     pkg, 
                     spec$version %||% "*", 
                     spec$source %||% "cran"))
        } else {
          cat(sprintf("  %s: %s\n", pkg, spec))
        }
      }
    }
  }
}

#' Sync dependencies: install from TOML and update lockfile
sync_dependencies <- function() {
  cat("Syncing dependencies from Rproject.toml...\n")
  
  # Install all dependency groups
  install_from_toml(c("dependencies", "dev-dependencies"))
  
  # Clean unused packages
  cat("Cleaning unused packages...\n")
  renv::clean(prompt = FALSE)
  
  cat("Sync complete!\n")
}

#' Null coalescing operator
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs) || length(lhs) == 0) rhs else lhs
}

# Helper function to convert R list to TOML string
toTOML <- function(data) {
  # This is a simplified TOML writer - in practice you'd want a more robust solution
  lines <- character()
  
  for (section in names(data)) {
    if (is.list(data[[section]])) {
      lines <- c(lines, paste0("[", section, "]"))
      
      for (key in names(data[[section]])) {
        value <- data[[section]][[key]]
        if (is.list(value)) {
          # Handle complex values like { source = "bioconductor" }
          parts <- sapply(names(value), function(k) {
            v <- value[[k]]
            if (is.character(v)) {
              paste0(k, ' = "', v, '"')
            } else {
              paste0(k, " = ", v)
            }
          })
          lines <- c(lines, paste0(key, " = { ", paste(parts, collapse = ", "), " }"))
        } else if (is.character(value)) {
          lines <- c(lines, paste0(key, ' = "', value, '"'))
        } else {
          lines <- c(lines, paste0(key, " = ", value))
        }
      }
      lines <- c(lines, "")
    }
  }
  
  return(lines)
}

# Command-line interface when run as script
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript renv_toml.R <command> [args]\n")
    cat("Commands:\n")
    cat("  install [groups]   - Install dependencies from TOML\n")
    cat("  add <pkg> [ver]    - Add new dependency\n") 
    cat("  remove <pkg>       - Remove dependency\n")
    cat("  list [group]       - List dependencies\n")
    cat("  sync               - Sync all dependencies\n")
    quit(status = 1)
  }
  
  command <- args[1]
  
  switch(command,
    "install" = {
      groups <- if (length(args) > 1) args[-1] else c("dependencies")
      install_from_toml(groups)
    },
    "add" = {
      if (length(args) < 2) stop("Package name required")
      package <- args[2]
      version <- if (length(args) > 2) args[3] else "*"
      add_dependency(package, version)
    },
    "remove" = {
      if (length(args) < 2) stop("Package name required")
      remove_dependency(args[2])
    },
    "list" = {
      group <- if (length(args) > 1) args[2] else NULL
      list_dependencies(group)
    },
    "sync" = {
      sync_dependencies()
    },
    stop("Unknown command: ", command)
  )
}