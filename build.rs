use chrono::Utc;
use git2::Repository;



fn main() {
    // Get the current time
    let now = Utc::now();
    println!("cargo:rustc-env=BUILD_TIME={}", now.to_rfc3339());

    // Get the current Git commit hash and compute its MD5 sum
    let repo_paths = [".", "../"];
    for path in repo_paths {
        if let Ok(repo) = Repository::open(path) {
            if let Ok(head) = repo.head() {
                if let Some(oid) = head.target() {
                    let commit_hash = oid.to_string();
                    let digest = md5::compute(commit_hash.clone());
                    let md5sum = format!("{:x}", digest);
                    println!("cargo:rustc-env=GIT_COMMIT_MD5={}", md5sum);
                    println!("cargo:warning=MD5 sum of the current Git commit: {}", md5sum);
                    return;
                }
            }
        }
    }

    println!("cargo:warning=Failed to open Git repository.");
}

